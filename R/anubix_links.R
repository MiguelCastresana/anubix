#' Create the Link Matrix Between Genes and Pathways
#'
#' @description
#' Computes the number or strength of connections that each gene in a network has to each pathway.
#' This matrix can be used for downstream network enrichment analysis.
#'
#' @usage
#' anubix_links(network, cores = 2, pathways, cutoff = 0.8, network_type = "weighted")
#'
#' @param network A data.frame representing the gene network.
#'   \itemize{
#'     \item If unweighted: Two columns (gene1, gene2), where each row represents an edge.
#'     \item If weighted: Three columns, with the third column containing edge weights.
#'   }
#'
#' @param cores Integer. Number of CPU cores to use for parallel computation. Default is 2.
#'
#' @param pathways A data.frame with two columns:
#'   \itemize{
#'     \item Column 1: gene names.
#'     \item Column 2: corresponding pathway names.
#'   }
#'
#' @param cutoff Numeric. Threshold to filter edges in a weighted network by confidence. Default is \strong{0.8}.
#'
#' @param network_type Character. Either \code{"weighted"} or \code{"unweighted"}. Default is \strong{"weighted"}.
#'
#' @return A numeric matrix where:
#' \itemize{
#'   \item Rows correspond to genes in the network.
#'   \item Columns correspond to pathways.
#'   \item Entries represent the number (or sum of weights) of links from each gene to each pathway.
#' }
#'
#' @import igraph
#' @importFrom purrr map set_names
#' @importFrom tibble as_tibble
#' @importFrom dplyr group_keys group_split pull
#' 
#' @export
#'
#' @seealso
#' \code{\link{anubix}},\code{\link{anubix_transitivity}}, \code{\link{example_anubix}}, \code{\link{anubix_clustering}}
#'
#' @examples
#' \dontrun{
#' # Example run on provided example data
#' anubix_links(
#'   network = example_anubix$network,
#'   pathways = example_anubix$pathway_set,
#'   cutoff = 0.8,
#'   network_type = "weighted"
#' )
#' }


anubix_links <- function(network, pathways, cutoff = 0.8, network_type = "weighted") {

  if (is.null(network)) {
    stop("Network file is missing", call. = FALSE)
  } else if (is.null(pathways)) {
    stop("Pathways file is missing", call. = FALSE)
  } else if (is.null(cutoff)) {
    cutoff <- 0.8
  }

  if (is.null(pathways) || ncol(pathways) != 2) {
    stop("Pathways missing or the file is not in a proper format.")
  }
  if (is.null(network) || ncol(network) < 2) {
    stop("A network is required or it is not in a proper format.")
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1) {
    stop("Link confidence cutoff is not in a proper format.")
  }
  if (is.null(network_type)) {
    network_type <- "weighted"
  }
  if (!network_type %in% c("weighted", "unweighted")) {
    stop("network_type must be either 'weighted' or 'unweighted'.")
  }
  if (network_type == "weighted" && ncol(network) < 3) {
    stop("A weighted network must contain a third column with edge weights.")
  }

  if (network_type == "weighted") {
    net <- network[network[, 3] >= cutoff, , drop = FALSE]
    net <- data.frame(
      from = net[, 1],
      to = net[, 2],
      weight = as.numeric(net[, 3])
    )
  } else {
    net <- data.frame(
      from = network[, 1],
      to = network[, 2]
    )
  }

  # convert network to igraph
  net_graph <- graph_from_data_frame(net, directed = FALSE, vertices = NULL)

  # get the adjacency matrix
  edge_attr <- if (network_type == "weighted") "weight" else NULL
  net_adjacency <- as_adjacency_matrix(net_graph, attr = edge_attr, names = TRUE)
  
  # add 2 in the diagonal, to boost overlap
  # overlap = diag(2,dim(net_adjacency)[1],dim(net_adjacency)[1])
  # net_adjacency = net_adjacency+ overlap
  # get all nodes in the graph
  net_nodes <- V(net_graph)$name

  # filter pathways for nodes not in FC
  pathways <- as_tibble(pathways)
  names(pathways) <- c("symbol", "pathway_names")
  pathway_list_filtered <- pathways %>%
    dplyr::filter(symbol %in% net_nodes)

  # convert pathway_df to a list
  pathway_names <-
    pathway_list_filtered %>% group_by(pathway_names) %>% group_keys() %>% pull(1)
  pathway_list_clean <-
    pathway_list_filtered %>% group_by(pathway_names) %>% group_split(.keep = FALSE) %>% set_names(pathway_names) %>% purrr::map(pull, symbol)
  
  # Create the links-matrix with zeros
  links_matrix <- matrix(0, nrow = length(net_nodes), ncol = length(pathway_list_clean),
                         dimnames = list(net_nodes, names(pathway_list_clean)))

  # Fill the matrix with the number or strength of links each gene has to a pathway.
  for (pathw in names(pathway_list_clean)) {
    path_nodes <- as.character(as.vector(pathway_list_clean[[pathw]]))
    links_matrix[, pathw] <- Matrix::rowSums(net_adjacency[, path_nodes, drop = FALSE])
  }
  links_matrix <- as.data.frame(links_matrix)
  return(links_matrix)
}

