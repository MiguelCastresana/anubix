
#' Perform ANUBIX transitivity test for the provided gene sets
#'
#' @note anubix_links() needs to be run beforehand.
#' @description Computes ANUBIX transitivity, an accurate test for network enrichment analysis between query sets and pathway sets. 
#' Instead of normal random sampling, it performs constrained random sampling, taking the degree of the nodes into account. 
#' Additionally, it incorporates the gene set’s transitivity to compute enrichment.
#' 
#' @usage anubix_transitivity(network, links_matrix, genesets, pathways, cores = 2, cutoff = 0.8,
#' sampling = 2000, network_type = "weighted")
#' 
#' @param network A data.frame. Two columns if the network has no weights, where the first column and second column are genes. Each row represents a link between genes. If weighted, the third column contains the weights of the links.
#' @param links_matrix A numeric matrix. Stores the links each gene has to each pathway. Rows correspond to genes in the network; columns correspond to pathways.
#' @param genesets A data.frame with two columns: the first column contains genes, the second column the experiment they belong to.
#' @param pathways A data.frame with two columns: the first column contains genes, the second column the pathway they belong to.
#' @param cores Numeric. Number of CPU cores used by the algorithm. Default is \strong{2}.
#' @param sampling Numeric. Number of random samplings to construct the null distribution. Default is \strong{2000}.
#' @param cutoff Numeric. Link confidence threshold for weights between genes. Default is \strong{0.8}.
#' @param network_type Character. Either "weighted" or "unweighted". Default is \strong{weighted}.
#' 
#' @importFrom TailRank dbb
#' @importFrom dplyr %>% select
#' @importFrom optimx optimr
#' @importFrom stats p.adjust
#' @importFrom igraph transitivity graph_from_edgelist
#' @import parallel
#' 
#' @export
#' 
#' @return A data.frame with columns:
#' \itemize{
#'   \item geneset   - Gene set under study.
#'   \item pathway   - Pathway under study.
#'   \item obv_links - Observed number of links between the gene set and the pathway.
#'   \item exp_mean  - Expected number of links between the gene set and the pathway.
#'   \item overlap   - Number of genes shared by the gene set and the pathway.
#'   \item p-value   - p-value of the test.
#'   \item q-value   - Corrected p-value using Benjamini-Hochberg procedure.
#'   \item FWER      - Corrected p-value using Bonferroni correction.
#' }
#' 
#' @seealso \code{\link{anubix_links}}, \code{\link{example_anubix}}, \code{\link{anubix_clustering}}, \code{\link{anubix}}
#' 
#' @examples
#' \dontrun{
#'  anubix_transitivity(
#'    network = example_anubix$network,
#'    links_matrix = example_anubix$links_genes,
#'    genesets = example_anubix$gene_set,
#'    pathways = example_anubix$pathway_set,
#'    cores = 2,
#'    cutoff = 0.8,
#'    sampling = 2000,
#'    network_type = "weighted"
#'  )
#' }



anubix_transitivity = function(network,links_matrix,
                              genesets,
                              pathways,
                              cores = 2, cutoff = 0.8,
                              sampling = 2000,network_type = "weighted",
                              callback = NULL,
                              website = FALSE) {
  if (is.null(network)){

    stop("Network file is missing", call.=FALSE)
  }
  else if (is.null(links_matrix)) {
    stop("Precomputed file for the links per gene is missing",
         call. = FALSE)
  }
  else if (is.null(genesets)) {
    stop("Query gene sets are missing", call. = FALSE)
  }
  else if (is.null(pathways)) {
    stop("Pathway file is missing", call. = FALSE)
  }
  else if (is.null(cores)) {
    cores = 2
  }
  else if(is.null(cutoff)){

    cutoff = 0.8
  }
  else if (is.null(sampling)) {
    sampling = 2000
  }
  if (is.null(links_matrix) | class(links_matrix) != "data.frame")
    stop("Please introduce the proper data for the link count for each gene.")
  if (class(cores) != "numeric" | detectCores() < cores)
    stop("Please introduce a proper value.")
  # genesets = as.data.frame(genesets)
  # if (ncol(genesets) < 2){
  #
  #   genesets[,2] = rep("geneset1",nrow(genesets))
  # }
  if (class(cutoff)!="numeric")
    stop("Link confidence cutoff is not in a proper format.")
  if (is.null(pathways) | class(pathways) != "data.frame" |
      ncol(pathways) != 2)
    stop("Pathways missing or the file is not in a proper format.")
  if (class(sampling) != "numeric")
    stop("Please introduce a correct value for the total number of random samplings.")
  if(is.null(network_type)){network_type = "weighted"}



  if(network_type=="weighted" ){
    net = network[which(network[,3]>=cutoff),]
  }else{

    net = network
  }


  # Transform data
  pathways = as.data.frame(pathways)
  pathways[,1] = as.vector(pathways[,1])
  pathways[,2] = as.vector(pathways[,2])

  genesets = as.data.frame(genesets)
  genesets[,1] = as.vector(genesets[,1])
  genesets[,2] = as.vector(genesets[,2])

  
  genesets = genesets[which(genesets[,1]%in%rownames(links_matrix)),]
  group_sets = unique(as.vector(genesets[, 2]))
  
  freq_genesets = as.data.frame(table(genesets[,2]))
  freq_genesets = freq_genesets[order(freq_genesets$Freq),]
  pos = which(freq_genesets[,2]<2)
  
  if (length(pos) > 0) {
    select_genesets <- as.vector(freq_genesets$Var1[-pos])
  } else {
    # If no groups have less than 2 elements, keep all groups
    select_genesets <- as.vector(freq_genesets$Var1)
  }
  
  genesets = genesets[which(genesets[,2]%in%select_genesets),]
  genesets[,2] = as.vector(genesets[,2])
  group_sets = unique(as.vector(genesets[, 2]))
  
  net = as.data.frame(net)


  pathways = pathways[which(pathways[, 2] %in% colnames(links_matrix)),
  ]
  group_sets = unique(as.vector(genesets[, 2]))
  group_paths = unique(as.vector(colnames(links_matrix)))
  length_genesets = numeric()
  i = 1
  for (i in 1:length(group_sets)) {
    sub = genesets[which(genesets[, 2] %in% group_sets[i]),
    ]
    length_genesets[i] = nrow(sub)
  }
  length_pathways = numeric()
  i = 1
  for (i in 1:length(group_paths)) {
    sub = pathways[which(pathways[, 2] %in% group_paths[i]),
    ]
    length_pathways[i] = nrow(sub)
  }
  degreeList = as.data.frame(table(c(as.vector(net[, 1]), as.vector(net[,
                                                                        2]))))
  degreeList_sorted = degreeList[with(degreeList, order(Freq)),
  ]
  d = unique(degreeList_sorted[, 2])
  i = 1
  deg_possibilities = list()
  for (i in 1:length(d)) {
    g = as.vector(degreeList_sorted[which((degreeList_sorted[,
                                                             2] == d[i])), 1])
    deg_possibilities[[i]] = g
  }
  names(deg_possibilities) = d
  l = sapply(deg_possibilities, function(x) length(x))
  p = 1
  c = 1
  i = 1
  pos = numeric()
  l_sum = 0
  times = 0
  new_deg_list = list()
  poss = list()
  names_deg = vector()
  for (i in 1:length(deg_possibilities)) {
    l_sum = l_sum + l[i]
    if (l_sum >= 100) {
      pos[p] = i
      p = p + 1
    }
    if (l_sum < 100) {
      times = times + 1
      pos[p] = i
      p = p + 1
      if (i == length(deg_possibilities)) {
        new_deg_list[[c]] = as.vector(unlist(deg_possibilities[pos[1]:i]))
        a = names(deg_possibilities[pos[1]])
        b = names(deg_possibilities[(i)])
        names_deg[c] = paste(a, b, sep = "_to_")
        poss[[c]] = as.numeric(names(deg_possibilities[pos[1]:i]))
      }
      next
    }
    new_deg_list[[c]] = as.vector(unlist(deg_possibilities[pos[1]:i]))
    a = names(deg_possibilities[pos[1]])
    b = names(deg_possibilities[(i)])
    names_deg[c] = paste(a, b, sep = "_to_")
    poss[[c]] = as.numeric(names(deg_possibilities[pos[1]:i]))
    c = c + 1
    times = 0
    l_sum = 0
    pos = numeric()
    p = 1
  }
  names(new_deg_list) = names_deg
  l_new = as.vector(sapply(new_deg_list, function(x) length(x)))
  l_pos = as.vector(sapply(poss, function(x) length(x)))
  i = 1
  sub1 = numeric()
  sub = numeric()
  for (i in 1:length(l_pos)) {
    sub = c(sub, poss[[i]])
    ok = rep(i, length(poss[[i]]))
    sub1 = c(sub1, ok)
  }
  map_pos_deg_list = as.data.frame(cbind(sub, sub1))
  "%!in%" <- function(x, y) !(x %in% y)
  sampling_generator = function(geneset, n_randomizations) {
    b = as.data.frame(table(as.character(as.vector(geneset[[2]]))))
    b[, 1] = as.numeric(as.vector(b[, 1]))
    b = b[with(b, order(Var1)), ]
    gens = as.vector(geneset[[1]])
    a = lapply(new_deg_list, function(x) x[which(x %!in%
                                                   gens)])
    pos = map_pos_deg_list[which(map_pos_deg_list[, 1] %in%
                                   as.character(b[, 1])), 2]
    times = c(1:n_randomizations)
    ok = as.data.frame(cbind(pos, b[, 2]))
    pp = unique(as.vector(ok[, 1]))
    newlist = list()
    i = 1
    for (i in 1:length(pp)) {
      ss1 = sum(ok[which(ok[, 1] %in% pp[i]), 2])
      ss2 = pp[i]
      newlist[[i]] = c(ss2, ss1)
    }
    sub = do.call(rbind.data.frame, newlist)
    r_genes = list()
    i = 1
    for (i in 1:sampling) {
      r_genes[[i]] = as.vector(unlist(mapply(FUN = function(x,
                                                            y) {
        set.seed(times[i])
        sample(x, y, replace = F)
      }, a[sub[, 1]], sub[, 2], SIMPLIFY = FALSE)))
    }
    return(r_genes)
  }

  links_geneset = function(geneset) {
    links = links_matrix[geneset, ]
    links = as.data.frame(links)
    links = sapply(links, function(x) sum(x))
    links = unname(links, force = FALSE)
    return(links)
  }
  index_genes = rownames(links_matrix)
  filtered_genesets = genesets[which(as.vector(genesets[, 1]) %in%
                                       as.vector(degreeList[, 1])), ]
  group_sets = unique(as.vector(filtered_genesets[, 2]))
  length_genesets_filtered = numeric()
  i = 1
  for (i in 1:length(group_sets)) {
    sub = filtered_genesets[which(filtered_genesets[, 2] %in%
                                    group_sets[i]), ]
    length_genesets_filtered[i] = nrow(sub)
  }
  
  if (length(length_genesets_filtered) > 5) {
    set = seq(from = 1, to = length(length_genesets_filtered), by = 5)
    set[length(set)] = length(length_genesets_filtered)
  }else {
    set = length(length_genesets)
  }
  
  i = 1
  deg = numeric()
  for (i in 1:nrow(filtered_genesets)) {
    deg[i] = as.numeric(as.vector(degreeList[which(degreeList[,
                                                              1] %in% as.vector(filtered_genesets[i, 1])), 2]))
  }
  filtered_genesets[, 3] = deg

  if(nrow(filtered_genesets)<1){
    result <- data.frame(matrix(ncol = 8, nrow = 0))
    names(result) = c("geneset",
                      "pathway",
                      "obv_links",
                      "exp_mean",
                      "overlap",
                      "p-value",
                      "q-value",
                      "FWER")

  }else{

    i = 1
    geneset_test_list = list()
    for (i in 1:length(group_sets)) {
      geneset_test = filtered_genesets[which(filtered_genesets[,
                                                               2] %in% group_sets[i]), c(1, 3)]
      geneset_test_list[[i]] = geneset_test
    }
    genesets = filtered_genesets
    groups = group_sets
    i = 1
    query_list = list()
    for (i in 1:length(groups)) {
      sub = filtered_genesets[which(filtered_genesets[,
                                                      2] %in% groups[i]), 1]
      pos = which(index_genes %in% sub)
      query_list[[i]] = pos
    }
    i = 1
    overlap_list = list()
    for (i in 1:length(groups)) {
      sub = as.vector(genesets[which(genesets[, 2] %in%
                                       groups[i]), 1])
      overlap_list[[i]] = sub
    }
    links_geneset_real = function(geneset) {
      links = links_matrix[geneset, ]
      links = as.data.frame(links)
      links = sapply(links, function(x) sum(x))
      links = unname(links, force = FALSE)
      return(links)
    }
    links_geneset = function(geneset) {
      links = links_matrix[which(rownames(links_matrix) %in%
                                   as.vector(unlist(geneset))), ]
      links = as.data.frame(links)
      links = sapply(links, function(x) sum(x))
      links = unname(links, force = FALSE)
      return(links)
    }



    loglik = function(inits, x) {
      A <- inits[1]
      B <- inits[2]
      Y = x[, 3]
      N = x[, 2]
      -sum(lgamma(abs(A) + abs(B)) - lgamma(abs(A)) -
             lgamma(abs(B)) + lgamma(Y + abs(A)) + lgamma(N -
                                                            Y + abs(B)) - lgamma(N + abs(A) + abs(B)))
    }
    stat_computation = function(c, observed1, max1) {
      n = max1
      m_1 = mean(sol[[c]])
      m_2 = mean(sol[[c]]^2)
      alpha = (n * m_1 - m_2)/(n * (m_2/m_1 - m_1 -
                                      1) + m_1)
      beta = (n - m_1) * (n - m_2/m_1)/(n * (m_2/m_1 -
                                               m_1 - 1) + m_1)
      inits = c(alpha, beta)
      if (any(is.nan(inits)) == TRUE) {
        pvalue = NA
      }
      else if (any(inits == 0) == TRUE) {
        pvalue = NA
      }
      else {
        dat = as.data.frame(cbind(1:times, rep(max1,
                                               times), sol[[c]]))
        optim.tas = optimx::optimr(par = inits, 
                                   fn = loglik,  x = dat, method = "L-BFGS-B", control = list(allmeth = "L-BFGS-B", 
                                                                                              allpkg = "stats"))
        optim.tas$par = abs(optim.tas$par)
        pvalue = 0.5 * dbb(observed1, max1, optim.tas$par[1],
                           optim.tas$par[2]) + sum(dbb((observed1 +
                                                          1):(max1), (max1), optim.tas$par[1], optim.tas$par[2]))
      }
      return(pvalue)
    }
    transitivity_f = function(x) {
      a = as.matrix(net[which(net[, 1] %in% x & net[,
                                                    2] %in% x), c(1, 2)])
      g <- graph_from_edgelist(a)
      mc <- transitivity(g, type = "undirected")
      tr = mean(mc, na.rm = TRUE)
      tr = tr
      return(tr)
    }


      total = 0
      timer = 1
      cc = 1
      result_t <- list()



      for (timer in 1:length(set)) {

        no_cores <- cores
        cl <- makeCluster(no_cores)
        times = sampling
        geneset_test_list1 = geneset_test_list[cc:set[timer]]
        data = parLapply(cl, geneset_test_list1, function(x) as.data.frame(x))
        clusterExport(cl, list("sampling_generator",
                               "%!in%", "new_deg_list", "map_pos_deg_list",
                               "times"), envir = environment())
        length_geneset1 = length_genesets_filtered[cc:set[timer]]
        group_sets1 = group_sets[cc:set[timer]]
        prueba1 = parSapply(cl, data, function(x, y) sampling_generator(x,
                                                                        y), y = times)
        prueba1 = lapply(prueba1, function(x) {
          as.vector(unlist(x))
        })
        clusterExport(cl, list("links_geneset", "links_matrix"),
                      envir = environment())
        query = parLapply(cl, prueba1, function(x) links_geneset(x))
        chunked_list <- split(prueba1, ceiling(seq_along(prueba1)/1))
        results_chance <- unlist(lapply(chunked_list,
                                        function(x) {
                                          lapply(x, function(y) {
                                            transitivity_f(y)
                                          })
                                        }))
        results_chance <- ifelse(is.nan(results_chance),
                                 0, results_chance)
        results_chance = results_chance + 1
        query <- lapply(seq_along(query), function(i) round((query[[i]] *
                                                               results_chance[i])))
        geneset_list <- split(genesets[, 1], genesets[,
                                                      2])
        geneset_list <- split(geneset_list, ceiling(seq_along(geneset_list)/1))
        results_real <- unlist(lapply(geneset_list, function(x) {
          lapply(x, function(y) {
            transitivity_f(y)
          })
        }))
        results_real <- ifelse(is.nan(results_real),
                               0, results_real)
        results_real = results_real + 1
        m <- length(query[[1]])
        clusterExport(cl, list("m", "query"), envir = environment())
        remove(prueba1)
        gc()
        information_list_true = parLapply(cl, 1:m, function(j) sapply(query,
                                                                      "[[", j))
        remove(query)
        t1 = Sys.time()
        real_genesets = lapply(query_list[cc:set[timer]],
                               function(x) links_geneset_real(x))
        real_genesets <- lapply(seq_along(real_genesets),
                                function(i) round((real_genesets[[i]] * results_real[i])))
        stopCluster(cl)

        no_cores <- cores
        cl <- makeCluster(no_cores)
        clusterExport(cl, list("information_list_true",
                               "real_genesets", "query_list", "links_matrix",
                               "times", "group_paths", "pathways", "length_pathways",
                               "length_geneset1"), envir = environment())
        sol = parSapply(cl, information_list_true, function(x) split(x,
                                                                     ceiling(seq_along(x)/times)))
        clusterExport(cl, list("sol"), envir = environment())
        remove(information_list_true)
        observed = as.vector(unlist(real_genesets))
        g_set = parLapply(cl, overlap_list[cc:set[timer]],
                          function(x) as.vector(unlist(x)))
        overlapp = parLapply(cl, g_set, function(x) as.vector(sapply(group_paths,
                                                                     function(y) length(x[which(x %in% as.vector(pathways[which(pathways[,
                                                                                                                                         2] %in% y), 1]))]))))
        overlapp = unlist(overlapp)
        remove(g_set)
        length_pathways = as.numeric(as.vector(length_pathways))
        length_genesets1 = as.numeric(as.vector(length_geneset1))
        max = unlist(parLapply(cl, length_geneset1, function(x) as.vector(sapply(length_pathways,
                                                                                 function(y) (y * x) - min(y, x)))))

        max =  2*(max - overlapp)
        stopCluster(cl)

        positions = seq(from = 1, to = length(length_pathways) *
                          length(length_geneset1), by = length(length_geneset1))
        positions = rep(positions, length(length_geneset1))
        summ = rep(0:(length(length_geneset1) - 1), each = length(length_pathways))
        positions = positions + summ
        n.cores = cores
        if (.Platform$OS.type == "windows") {
          n.cores = 1
        }
        pvalues = as.vector(unlist(mcmapply(function(x,
                                                     y, z) stat_computation(x, y, z), positions,
                                            observed, max, mc.cores = n.cores)))

        expected = as.numeric(unlist(lapply(sol, function(x) mean(x))))
        remove(sol)
        expected = expected[positions]
        result = as.data.frame(cbind(rep(group_sets1,
                                         each = length(length_pathways)), rep(group_paths,
                                                                              length(length_genesets1)), observed, expected,
                                     overlapp, pvalues))
        result_t[[timer]] = result
        remove(result)
        remove(pvalues)
        remove(expected)
        remove(positions)
        remove(overlapp)
        cc = set[timer] + 1
        print(timer)
        if (!is.null(callback)) {
          callback("100% of the process done")
        }

        # if(timer>1){break}
      }
      result_t = do.call(rbind, result_t)
      result_t[, 3] = as.numeric(as.vector(result_t[, 3]))
      result_t[, 4] = as.numeric(as.vector(result_t[, 4]))
      result_t[, 5] = as.numeric(as.vector(result_t[, 5]))
      result_t[, 6] = as.numeric(as.vector(result_t[, 6]))
      result_t = result_t[complete.cases(result_t),]
      result_t[, 7] = p.adjust(result_t[, 6], method = "BH")
      result_t[, 8] = p.adjust(result_t[, 6], method = "bonferroni")
      names(result_t) = c("geneset",
                          "pathway",
                          "obv_links",
                          "exp_mean",
                          "overlap",
                          "p-value",
                          "q-value",
                          "FWER")
    

    return(result_t)

  }



}
