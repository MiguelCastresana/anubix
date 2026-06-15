script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
local_library <- file.path(repo_root, ".Rlib")

if (dir.exists(local_library)) {
  .libPaths(c(local_library, .libPaths()))
}

setwd(repo_root)
pkgload::load_all(repo_root, export_all = FALSE, helpers = FALSE)

cache_dir <- file.path(repo_root, "tests", "testthat", "_cache")
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
selected_pathways <- names(pathway_sizes[pathway_sizes >= 8])[1:8]
selected_pathways <- selected_pathways[!is.na(selected_pathways)]
pathways <- pathways[pathways$pathway %in% selected_pathways, , drop = FALSE]

requested_source <- Sys.getenv("ANUBIX_SOURCE_PATHWAY", unset = "")
source_pathway <- if (nzchar(requested_source)) requested_source else selected_pathways[[1]]
if (!source_pathway %in% selected_pathways) {
  selected_pathways <- unique(c(source_pathway, selected_pathways))
  pathways <- pathways[pathways$pathway %in% selected_pathways, , drop = FALSE]
}
query_genes <- unique(pathways$gene[pathways$pathway == source_pathway])
query_genes <- query_genes[query_genes %in% nodes]
query_genes <- head(query_genes, 12)
genesets <- data.frame(gene = query_genes, geneset = paste0("query_from_", source_pathway))

links_matrix <- anubix_links(
  network = network,
  pathways = pathways,
  network_type = "unweighted"
)

set.seed(1)
result <- anubix(
  network = network,
  links_matrix = links_matrix,
  genesets = genesets,
  pathways = pathways,
  cores = 1,
  sampling = 200,
  network_type = "unweighted"
)

result <- result[order(result[["p-value"]], result[["q-value"]]), ]
result$source_pathway <- result$pathway == source_pathway

print(list(
  network_edges = nrow(network),
  network_nodes = length(nodes),
  pathway_rows = nrow(pathways),
  source_pathway = source_pathway,
  query_size = nrow(genesets),
  selected_pathways = selected_pathways
))

print(result[, c(
  "geneset",
  "pathway",
  "obv_links",
  "exp_mean",
  "overlap",
  "p-value",
  "q-value",
  "FWER",
  "source_pathway"
)], row.names = FALSE)
