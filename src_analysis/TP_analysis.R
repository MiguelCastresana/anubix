# Load necessary libraries
library(parallel)
library(fastmatch)
library(TailRank)
library(optimr)

# Helper to get a list of genes from KEGGA split by pathway
prepare_gene_sets <- function(file_path) {
  df <- read.delim(file_path, header = TRUE)
  split(df, f = df[, 2])
}

# Compute sum of links per geneset in link matrix
links_geneset_real <- function(geneset, output_matrix) {
  subset <- output_matrix[rownames(output_matrix) %fin% as.vector(unlist(geneset[1])), ]
  links <- colSums(subset)
  unname(links)
}

# Function to generate random genes
get_random_string <- function(n = 1, length = 15) {
  replicate(n, paste(sample(c(0:9, letters, LETTERS), length, replace = TRUE), collapse = ""))
}

# Sample genes from total universe
samples <- function(dat_length, genesall) {
  sample(genesall, dat_length, replace = FALSE)
}

# Beta-binomial log-likelihood
loglik <- function(inits, dat) {
  A <- inits[1]
  B <- inits[2]
  Y <- dat[, 3]
  N <- dat[, 2]
  -sum(
    lgamma(abs(A) + abs(B)) - lgamma(abs(A)) - lgamma(abs(B)) +
      lgamma(Y + abs(A)) + lgamma(N - Y + abs(B)) - lgamma(N + abs(A) + abs(B))
  )
}

# Run ANUBIX evaluation
run_anubix_benchmark <- function(kegga_path, keggb_path, link_matrix_path, net_path, output_path, times = 2000) {
  KEGGA <- read.delim(kegga_path, header = TRUE)
  KEGGB <- read.delim(keggb_path, header = TRUE)
  output1 <- read.delim(link_matrix_path, header = TRUE)
  net <- read.delim(net_path, header = TRUE)
  
  paths <- unique(as.vector(KEGGB[, 2]))
  geneset_test_list <- lapply(paths, function(p) KEGGA[which(KEGGA[, 2] %fin% p), 1])
  
  KGG_A <- prepare_gene_sets(kegga_path)
  KGG_B <- prepare_gene_sets(keggb_path)
  
  length_genesetA <- sapply(KGG_A, nrow)
  length_genesetB <- sapply(KGG_B, nrow)
  
  real_genesets <- lapply(KGG_A, function(x) links_geneset_real(x, output1))
  
  genes <- unique(c(as.vector(net[, 1]), as.vector(net[, 2])))
  genesall <- c(genes, get_random_string(20000 - length(genes), 15))
  
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores)
  clusterExport(cl, list("samples", "genesall", "times"), envir = environment())
  
  data_lengths <- sapply(geneset_test_list, length)
  prueba1 <- parSapply(cl, data_lengths, function(x) lapply(1:times, function(y) samples(x, genesall)))
  prueba1 <- lapply(prueba1, unlist)
  
  clusterExport(cl, list("links_geneset", "output1"), envir = environment())
  links_geneset <- function(geneset) {
    subset <- output1[rownames(output1) %fin% geneset, ]
    colSums(subset)
  }
  
  query <- parLapply(cl, prueba1, function(x) links_geneset(x))
  m <- length(query[[1]])
  clusterExport(cl, list("m", "query"), envir = environment())
  information_list_true <- parLapply(cl, 1:m, function(j) sapply(query, "[[", j))
  
  stopCluster(cl)
  
  pvalue <- vector()
  exp_mean <- vector()
  obv_links <- vector()
  geneset_c <- vector()
  geneset2 <- vector()
  y <- 1
  d <- times
  k <- 1
  
  for (j in seq_along(real_genesets)) {
    obv <- as.numeric(as.vector(real_genesets[[j]][j]))
    subset <- information_list_true[[j]][y:d]
    g_set <- as.vector(unlist(geneset_test_list[[j]][1]))
    overlap <- length(g_set[g_set %in% KEGGB[KEGGB[, 2] %in% paths[j], 1]])
    max_possible <- (length_genesetB[j] * length_genesetA[j]) - overlap
    
    y <- y + times
    d <- d + times
    
    if (mean(subset) == 0) subset[length(subset)] <- 1
    
    m_1 <- mean(subset)
    m_2 <- mean(subset^2)
    
    alpha <- abs((max_possible * m_1 - m_2) / (max_possible * (m_2 / m_1 - m_1 - 1) + m_1))
    beta <- abs((max_possible - m_1) * (max_possible - m_2 / m_1) / (max_possible * (m_2 / m_1 - m_1 - 1) + m_1))
    
    if (any(c(alpha, beta) <= 0)) next
    
    dat <- data.frame(idx = 1:times, N = rep(max_possible, times), Y = subset)
    fit <- tryCatch(optimr(c(alpha, beta), loglik, dat = dat), error = function(e) NULL)
    if (is.null(fit) || any(is.na(fit$par)) || fit$value < 0) next
    
    fit$par <- abs(fit$par)
    pval <- 0.5 * dbb(obv, max_possible, fit$par[1], fit$par[2]) +
      sum(dbb((obv + 1):(max_possible), max_possible, fit$par[1], fit$par[2]))
    
    pvalue[k] <- pval
    exp_mean[k] <- m_1
    obv_links[k] <- obv
    geneset_c[k] <- names(real_genesets[j])
    geneset2[k] <- paths[j]
    k <- k + 1
  }
  
  A <- data.frame(geneset = geneset_c, pathway = geneset2, obv_links = obv_links,
                  exp_mean = exp_mean, p_value = pvalue)
  A <- A[complete.cases(A), ]
  A$q_value <- p.adjust(A$p_value, method = "BH")
  
  write.table(A, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  return(A)
}