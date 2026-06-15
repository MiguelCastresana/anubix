test_that("ANUBIX runs with downloaded FunCoup and KEGG data", {
  skip_if_not(
    identical(Sys.getenv("ANUBIX_RUN_INTEGRATION"), "true"),
    "Set ANUBIX_RUN_INTEGRATION=true to run external-data integration tests."
  )

  cache_dir <- testthat::test_path("_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  network_file <- file.path(cache_dir, "FC6.0_M.jannaschii_compact.gz")
  if (!file.exists(network_file)) {
    download.file(
      "https://funcoup.org/download/network&FC6.0_M.jannaschii_compact.gz",
      network_file,
      mode = "wb",
      quiet = TRUE
    )
  }

  network_raw <- read.delim(gzfile(network_file), check.names = FALSE)
  network <- data.frame(
    gene_a = network_raw[["0:ProteinA"]],
    gene_b = network_raw[["1:ProteinB"]],
    weight = network_raw[["5:PPV"]]
  )
  network <- network[network$weight >= 0.90, , drop = FALSE]
  nodes <- unique(c(network$gene_a, network$gene_b))

  pathway_links <- read.delim(
    "http://rest.kegg.jp/link/pathway/mja",
    header = FALSE,
    col.names = c("kegg_gene", "pathway")
  )
  uniprot_links <- read.delim(
    "http://rest.kegg.jp/conv/uniprot/mja",
    header = FALSE,
    col.names = c("kegg_gene", "uniprot")
  )

  uniprot_links$uniprot <- sub("^up:", "", uniprot_links$uniprot)
  pathway_links$pathway <- sub("^path:", "", pathway_links$pathway)

  pathways <- merge(pathway_links, uniprot_links, by = "kegg_gene")
  pathways <- unique(pathways[, c("uniprot", "pathway")])
  names(pathways) <- c("gene", "pathway")
  pathways <- pathways[pathways$gene %in% nodes, , drop = FALSE]

  pathway_sizes <- sort(table(pathways$pathway), decreasing = TRUE)
  selected_pathways <- names(pathway_sizes[pathway_sizes >= 5])[1:4]
  selected_pathways <- selected_pathways[!is.na(selected_pathways)]
  pathways <- pathways[pathways$pathway %in% selected_pathways, , drop = FALSE]

  expect_gte(nrow(network), 100)
  expect_gte(length(selected_pathways), 2)
  expect_gte(nrow(pathways), 10)

  links_matrix <- anubix_links(
    network = network,
    pathways = pathways,
    network_type = "unweighted"
  )

  expect_s3_class(links_matrix, "data.frame")
  expect_equal(colnames(links_matrix), selected_pathways)
  expect_true(all(rowSums(links_matrix) >= 0))

  query_genes <- unique(pathways$gene[pathways$pathway == selected_pathways[[1]]])
  query_genes <- query_genes[query_genes %in% rownames(links_matrix)]
  query_genes <- head(query_genes, 3)
  genesets <- data.frame(gene = query_genes, geneset = "query_from_kegg")

  result <- anubix(
    network = network,
    links_matrix = links_matrix,
    genesets = genesets,
    pathways = pathways,
    cores = 1,
    sampling = 5,
    network_type = "unweighted"
  )

  expect_s3_class(result, "data.frame")
  expect_named(
    result,
    c("geneset", "pathway", "obv_links", "exp_mean", "overlap", "p-value", "q-value", "FWER")
  )
})
