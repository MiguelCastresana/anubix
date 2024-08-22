#' Creates the link_matrix data
#'
#'
#' @description Computation of all the links that each gene in the network has to each of the pathways.
#' @usage anubix_links(network,cores = 2,pathways,cutoff = 0.75,network_type = "weighted")
#' @param network A data.frame. Two columns if the network has no weights. Where column 1 and column 2 are genes. Each row means a link between genes. If the network is weighted, then the third column are the weights of the links.
#' @param pathways A data.frame of two columns. Column 1 are the genes and second column the pathway where they belong to.
#' @param cutoff A numeric value. Cutoff is defined as the link confidence threshold of the weights between genes. The default value is \strong{0.75}.
#' @param network_type Either "weighted" or "unweighted". The default value is \strong{weighted}.
#' @import igraph
#' @importFrom purrr map
#' @importFrom purrr set_names
#' @importFrom tibble as_tibble
#' @importFrom dplyr group_keys
#' @importFrom dplyr group_split
#' @importFrom dplyr pull
#' @export
#' @return A matrix with the following dimensions:
#' \itemize{
#'   \item nrows   -   Total number of genes in the studied network.
#'   \item ncols   -   Total number of pathways under study.
#' }
#'
#'
#' @seealso \code{\link{anubix}},\code{\link{example_anubix}},\code{\link{anubix_clustering}}
#'
#' @examples
#'
#' # Example with a tiny network:
#'\dontrun{
#' anubix_links(example_anubix$network,example_anubix$pathway_set,cutoff = 0.75, "weighted")
#' }


anubix_links = function(network,pathways,cutoff = 0.75,network_type = "weighted"){

  if (is.null(network)){

    stop("Network file is missing", call.=FALSE)
  }else if(is.null(pathways)){

    stop("Pathways file is missing", call.=FALSE)
  }else if(is.null(cutoff)){

    cutoff = 0.75
  }


  if (is.null(pathways) | ncol(pathways) != 2)
    stop("Pathways missing or the file is not in a proper format.")
  if (is.null(network) | ncol(network) < 2)
    stop("A network is required or it is not in a proper format.")
  if (class(cutoff)!="numeric")
    stop("Link confidence cutoff is not in a proper format.")
  if(is.null(network_type)){network_type = "weighted"}



  if(network_type=="weighted" ){
    net = network[which(network[,3]>=cutoff),]
  }else{

    net = network
  }
  # convert network to igraph
  net_graph <- graph_from_data_frame(net, directed = FALSE, vertices = NULL)
  # get the adjacency matrix
  net_adjacency <- as_adjacency_matrix(net_graph, names = TRUE)
  
  # add 2 in the diagonal, to boost overlap
  # overlap = diag(2,dim(net_adjacency)[1],dim(net_adjacency)[1])
  # net_adjacency = net_adjacency+ overlap
  # get all nodes in the graph
  net_nodes <- V(net_graph)$name

  # filter pathways for nodes not in FC
  pathways = as_tibble(pathways)
  names(pathways) = c("symbol","pathway_names")
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

  # Fill the start_matrix with number of links each gene has to a pathway
  for (pathw in names(pathway_list_clean)) {
    path_nodes <- as.character(as.vector(pathway_list_clean[[pathw]]))
    links_matrix[, pathw] <- Matrix::rowSums(net_adjacency[, path_nodes,drop=F])
  }
  links_matrix = as.data.frame(links_matrix)
  return(links_matrix)


}







