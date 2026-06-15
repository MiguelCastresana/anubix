# Benchmark MSigDB geneset against KEGG pathways using ANUBIX with convergence check

source(file.path("src_analysis", "paths.R"))

# Libraries
library(fastmatch)
library(parallel)
library(stringi)
library(TailRank)
library(optimr)

# --- 1. Load data ---
msigdb <- read.delim(benchmark_file("data", "msigdb"), header = TRUE)
geneset1 <- subset(msigdb, V2 == "DAIRKEE_CANCER_PRONE_RESPONSE_BPA")
geneset1$V2 <- as.character(geneset1$V2)

net <- read.delim(benchmark_file("data", "fc3.tsv"), header = TRUE)
genes <- unique(c(net[,1], net[,2]))

# --- 2. Build gene universe ---
get_random_string <- function(n = 1, length = 15) {
  replicate(n, paste0(sample(c(0:9, letters, LETTERS), length, replace = TRUE), collapse = ""))
}
rnd <- get_random_string(20000 - length(genes))
genesall <- c(genes, rnd)

output1 <- read.delim(benchmark_file("data", "link_matrix.tsv"), row.names = 1, header = TRUE)

# --- 3. Prepare KEGG pathways ---
KEGG <- read.delim(benchmark_file("data", "KEGG_pathways"), header = FALSE)
paths <- unique(KEGG$V2)
KGG_list <- split(KEGG, f = KEGG$V2)
length_genesetB <- sapply(KGG_list, nrow)

# Real link counts for each pathway
links_geneset_real <- function(geneset) {
  mat <- output1[rownames(output1) %fin% geneset$V1, ]
  colSums(mat)
}
real_genesets <- lapply(KGG_list, links_geneset_real)
length_geneset <- nrow(geneset1)

# Sampling function
samples_fun <- function(n) sample(genesall, n, replace = FALSE)

# --- 4. Empirical p-value sampling ---
times <- 2000
reps <- 100
pvalue_list <- vector("list", reps)

for (r in seq_len(reps)) {
  # parallel cluster
  cores <- max(1, detectCores() - 1)
  cl <- makeCluster(cores)
  clusterEvalQ(cl, { library(stringi); library(fastmatch) })
  clusterExport(cl, c("samples_fun", "genesall", "times"), envir = environment())
  
  # generate random gene sets
  rand_sets <- parSapply(cl, length_geneset, function(n) replicate(times, samples_fun(n), simplify = FALSE))
  rand_sets <- lapply(rand_sets, unlist)
  
  # count links in each random set
  links_rand <- function(gs) {
    mat <- output1[rownames(output1) %fin% gs, ]
    colSums(mat)
  }
  clusterExport(cl, c("output1", "links_rand"), envir = environment())
  query <- parLapply(cl, rand_sets, links_rand)
  
  # transpose to distributions per link
  m <- length(query[[1]])
  clusterExport(cl, c("query", "m"), envir = environment())
  info_list <- parLapply(cl, seq_len(m), function(j) sapply(query, `[[`, j))
  stopCluster(cl)
  
  # fit beta-binomial and compute p-values
  loglik <- function(params, dat) {
    A <- abs(params[1]); B <- abs(params[2])
    Y <- dat$Y; N <- dat$N
    -sum(lgamma(A+B) - lgamma(A) - lgamma(B) + lgamma(Y+A) + lgamma(N - Y + B) - lgamma(N + A + B))
  }
  
  s_vec <- seq(1, s <- times)
  y <- 1; d <- times; idx <- 1
  pval <- numeric()
  
  for (i in seq_along(real_genesets)) {
    for (j in seq_along(info_list)) {
      if (j > length(real_genesets[[i]])) break
      obs <- real_genesets[[i]][j]
      subset_vals <- info_list[[j]][y:d]
      
      overlap <- sum(geneset1$V1 %in% KEGG$V1[KEGG$V2 == paths[j]])
      max_links <- (length_genesetB[j] * length_geneset) - overlap
      
      if (mean(subset_vals) == 0) subset_vals[s] <- 1
      m1 <- mean(subset_vals); m2 <- mean(subset_vals^2)
      alpha <- abs((max_links * m1 - m2) / (max_links * (m2/m1 - m1 - 1) + m1))
      beta <- abs((max_links - m1) * (max_links - m2/m1) / (max_links * (m2/m1 - m1 - 1) + m1))
      if (alpha <=0 || beta <= 0) next
      
      dat <- data.frame(N = rep(max_links, times), Y = subset_vals)
      fit <- tryCatch(optimr(c(alpha, beta), loglik, dat = dat), error = function(e) NULL)
      if (is.null(fit) || fit$value < 0) next
      
      a <- abs(fit$par[1]); b <- abs(fit$par[2])
      ptemp <- 0.5 * dbb(obs, max_links, a, b) + sum(dbb((obs+1):max_links, max_links, a, b))
      pval[idx] <- ptemp
      idx <- idx + 1
    }
    y <- y + s; d <- d + s
  }
  pvalue_list[[r]] <- pval
}

# --- 5. Confidence intervals ---
m <- length(pvalue_list[[1]])
p_each <- lapply(seq_len(m), function(j) sapply(pvalue_list, `[[`, j))

normConfInt <- function(x, alpha=0.05) {
  m <- mean(x); sdv <- sd(x); n <- length(x)
  m + qt(1 - alpha/2, n-1) * sdv / sqrt(n) * c(-1, 1)
}
conf_int95 <- t(sapply(p_each, normConfInt))
conf_int95 <- data.frame(pathway = paths[seq_len(nrow(conf_int95))],
                         Lower_bound = conf_int95[,1], Upper_bound = conf_int95[,2])

# --- 6. Convergence: CV <= 0.02 ---
no_cores <- max(1, detectCores() - 1)
cl <- makeCluster(no_cores)
clusterEvalQ(cl, { library(stringi); library(fastmatch) })
clusterExport(cl, c("samples_fun","genesall"), envir = environment())

times_cv <- 12000
rand_sets_cv <- parSapply(cl, length_geneset, function(n) replicate(times_cv, samples_fun(n), simplify=FALSE))
rand_sets_cv <- lapply(rand_sets_cv, unlist)
clusterExport(cl, c("links_rand","output1"), envir = environment())
query_cv <- parLapply(cl, rand_sets_cv, links_rand)
m_cv <- length(query_cv[[1]])
clusterExport(cl, c("query_cv","m_cv"), envir = environment())
info_null <- parLapply(cl, seq_len(m_cv), function(j) sapply(query_cv, `[[`, j))
stopCluster(cl)

cv_cutoff <- 0.02
min_samples <- integer(length(info_null))
for (j in seq_along(info_null)) {
  cv <- 1; n <- 400
  while (cv > cv_cutoff && n <= times_cv) {
    sampm <- replicate(100, sample(info_null[[j]], n, FALSE))
    mvals <- colMeans(sampm)
    cv <- sd(mvals)/mean(mvals)
    if (mean(mvals) < 1) break
    n <- n + 100
  }
  min_samples[j] <- n
}

# Return outputs
data.frame(conf_int95, min_samples = min_samples)
