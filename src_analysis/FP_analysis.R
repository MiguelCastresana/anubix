# False positive benchmarking of gene sets in FunCoup network using ANUBIX null model

source(file.path("src_analysis", "paths.R"))

# --- 1. Libraries ---
library(fastmatch)
library(parallel)
library(stringi)
library(magrittr)
library(TailRank)
library(optimr)

# --- 2. Load network and prepare gene universe ---
net_full <- read.delim(benchmark_file("data", "fc3.tsv"), header = TRUE)
net <- subset(net_full, V3 >= 0.75)[, 1:3]
genes <- unique(c(net_full$V1, net_full$V2))

# Extend universe with random strings
get_random_string <- function(n, length = 12) {
  replicate(n, paste0(sample(c(0:9, letters, LETTERS), length, replace = TRUE), collapse = ""))
}
rnd <- get_random_string(20000 - length(genes))
genesall <- c(genes, rnd)

# Precompute node degrees in thresholded network
deg_df <- as.data.frame(table(c(net$V1, net$V2)), stringsAsFactors = FALSE)
names(deg_df) <- c("gene", "degree")

# --- 3. Load gene sets and pathways ---
# Random gene sets for false positive testing
KEGGA <- read.delim(benchmark_file("data", "randomsets_10000"), header = TRUE)
# KEGG pathways
KEGGB <- read.delim(benchmark_file("data", "KEGG_pathways"), header = FALSE)

sets <- unique(KEGGA$V2)
paths <- unique(KEGGB$V2)

# Build test list: one vector of genes per set
geneset_test_list <- lapply(sets, function(s) KEGGA$V1[KEGGA$V2 == s])
names(geneset_test_list) <- sets

# Real link counts function
tally_links <- function(gs, link_mat) {
  subm <- link_mat[rownames(link_mat) %fin% gs, , drop = FALSE]
  colSums(subm)
}

# Load ANUBIX link matrix
link_mat <- read.delim(benchmark_file("data", "link_matrix.tsv"), row.names = 1, header = TRUE)

# Prepare real links for each set
real_links_list <- lapply(geneset_test_list, tally_links, link_mat)
length_A <- lengths(geneset_test_list)
length_B <- table(KEGGB$V2)[paths]

# --- 4. Load precomputed null model samples ---
# Assumes `false_positive_test_sets_for_nullmodel` provides `prueba1`: list of random gene vectors
load(benchmark_file("data", "false_positive_test_sets_for_nullmodel"))
# `prueba1` should be a list of length equal to unique set size (110), each containing `times` random samples
prueba1 <- lapply(prueba1, unlist)

# --- 5. Build null distributions of link counts ---
links_rand <- function(gs) {
  subm <- link_mat[rownames(link_mat) %fin% gs, , drop = FALSE]
  colSums(subm)
}
no_cores <- max(1, detectCores() - 2)
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(stringi)); clusterEvalQ(cl, library(fastmatch))
clusterExport(cl, c("links_rand", "link_mat"), envir = environment())
query_null <- parLapply(cl, prueba1, links_rand)
stopCluster(cl)

# Transpose null distributions: one vector per pathway index
m <- length(query_null[[1]])
null_dists <- lapply(seq_len(m), function(j) sapply(query_null, `[[`, j))

# --- 6. Fit beta-binomial and compute p-values per pathway ---
loglik <- function(params, dat) {
  A <- abs(params[1]); B <- abs(params[2])
  Y <- dat$Y; N <- dat$N
  -sum(lgamma(A+B) - lgamma(A) - lgamma(B) +
         lgamma(Y+A) + lgamma(N-Y+B) - lgamma(N+A+B))
}

times <- length(prueba1[[1]])
pvalue_list <- vector("list", length(paths))

for (pi in seq_along(paths)) {
  pv <- numeric(length(sets))
  for (si in seq_along(sets)) {
    obs <- real_links_list[[si]][pi]
    subset_vals <- null_dists[[pi]]
    max_links <- length_A[si] * length_B[pi] -
      length(intersect(geneset_test_list[[si]], KEGGB$V1[KEGGB$V2 == paths[pi]]))
    if (mean(subset_vals) == 0) subset_vals[times] <- 1
    m1 <- mean(subset_vals); m2 <- mean(subset_vals^2)
    alpha <- abs((max_links*m1 - m2) / (max_links*(m2/m1 - m1 - 1) + m1))
    beta  <- abs((max_links-m1) * (max_links - m2/m1) / (max_links*(m2/m1 - m1 - 1) + m1))
    if (alpha <=0 || beta<=0) next
    dat <- data.frame(N = rep(max_links, times), Y = subset_vals)
    fit <- tryCatch(optimr(c(alpha,beta), loglik, dat=dat), error = function(e) NULL)
    if (is.null(fit) || any(is.na(fit$par)) || fit$value < 0) next
    a <- abs(fit$par[1]); b <- abs(fit$par[2])
    pv[si] <- 0.5*dbb(obs, max_links, a, b) + sum(dbb((obs+1):max_links, max_links, a, b))
  }
  pvalue_list[[pi]] <- pv
  message("Completed pathway ", pi, " of ", length(paths))
}
names(pvalue_list) <- paths

# The `pvalue_list` contains vectors of p-values for each false-positive set across pathways.

# Return or save as needed
pvalue_list
