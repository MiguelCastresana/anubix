
#' Perform ANUBIX for the provided gene sets
#'
#' @note anubix_links() needs to be run beforehand
#' @description Computes ANUBIX, an accurate test for network enrichment analysis between query sets and pathway sets. Instead of normal random sampling it does constrained random sampling, taking the degree of the nodes into account.
#' @usage anubix(network,links_matrix, genesets, pathways, cores = 2,cutoff = 0.8,
#' sampling = 2000,network_type = "weighted")
#' @param network A data.frame. Two columns if the network has no weights. Where column 1 and column 2 are genes. Each row means a link between genes. If the network is weighted, then the third column are the weights of the links.
#' @param links_matrix A numeric matrix. Links_matrix stores the links each gene has to each pathway. The number of rows equals to the number of genes in the network and the number of columns equals to the number of pathways under study.
#' @param genesets A data.frame of two columns. Column 1, are the genes and column 2 the experiment where they belong to.
#' @param pathways A data.frame of two columns. Column 1 are the genes and second column the pathway where they belong to.
#' @param cores A numeric value. Cores is defined as the number of cores used by the algorithm. The default value is \strong{2}
#' @param sampling A numeric value. Sampling is defined as the number of random samplings required to construct the null distribution. The default value is \strong{2000}
#' @param A numeric value. Cutoff is defined as the link confidence threshold of the weights between genes. The default value is \strong{0.8}.
#' @param network_type Either "weighted" or "unweighted". The default value is \strong{weighted}.
#' @importFrom TailRank dbb
#' @importFrom dplyr %>% select
#' @importFrom optimr optimr
#' @importFrom stats p.adjust
#' @import parallel
#' @export
#' @return A data frame with the following columns:
#' \itemize{
#'   \item geneset   -   Gene set under study.
#'   \item pathway   -   Pathway under study.
#'   \item obv_links -   Observed number of links between the gene set and the pathway.
#'   \item exp_mean  -   Expected number of links between the gene set and the pathway.
#'   \item overlap   -   Number of genes shared by the gene set and the pathway
#'   \item p-value   -   p-value of the test.
#'   \item q-value   -   Corrected p-value using Benjamini-Hochberg.
#'   \item FWER      -   Corrected p-value using Bonferroni correction.
#' }
#'
#' @seealso \code{\link{anubix_links}},\code{\link{anubix_transitivity}},\code{\link{anubix_clustering}},\code{\link{example_anubix}}
#'
#'
#'
#' @examples
#' # We provide with example data to be able to run ANUBIX.
#' \dontrun{
#'  anubix(network = example_anubix$network,links_matrix = example_anubix$links_genes,genesets = example_anubix$gene_set,
#'  pathways = example_anubix$pathway_set,cores = 2, cutoff = 0.8, sampling = 2000,network_type = "weighted")
#' }



anubix = function(network,links_matrix,
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

  ###
  pathways = pathways[which(pathways[,2]%in%colnames(links_matrix)),]
  group_sets = unique(as.vector(genesets[, 2]))
  group_paths = unique(as.vector(colnames(links_matrix)))
  length_genesets = numeric()
  i = 1
  for (i in 1:length(group_sets)) {
    sub = genesets[which(genesets[, 2] %in% group_sets[i]),]
    length_genesets[i] = nrow(sub)
  }
  length_pathways = numeric()
  i = 1
  for (i in 1:length(group_paths)) {
    sub = pathways[which(pathways[, 2] %in% group_paths[i]),]
    length_pathways[i] = nrow(sub)
  }



  degreeList = as.data.frame(table(c(as.vector(net[,1]),as.vector(net[,2]))))

  degreeList_sorted = degreeList[with(degreeList, order(Freq)), ]
  # through all the degrees in the net
  d = unique(degreeList_sorted[,2])
  i = 1
  deg_possibilities = list()
  for(i in 1:length(d)){


    g = as.vector(degreeList_sorted[which((degreeList_sorted[,2] == d[i])),1])

    deg_possibilities[[i]] = g

  }
  names(deg_possibilities) = d


  l = sapply(deg_possibilities, function(x) length(x))


  # We are going to make bins of at least 100 genes, to not to pick the same genes all the time



  p = 1
  c = 1
  i = 1
  pos = numeric()
  l_sum = 0
  times = 0
  new_deg_list = list()
  poss = list()
  names_deg = vector()
  for(i in 1:length(deg_possibilities)){


    l_sum = l_sum + l[i]
    if(l_sum>=100){
      pos[p] = i
      p = p + 1
    }
    if(l_sum<100){
      times = times + 1
      pos[p] = i
      p = p + 1
      if(i == length(deg_possibilities)){

        new_deg_list[[c]] = as.vector(unlist(deg_possibilities[pos[1]:i]))
        a = names(deg_possibilities[pos[1]]) # first degree
        b = names(deg_possibilities[(i)]) # last degree range
        names_deg[c] = paste(a,b,sep = "_to_")
        poss[[c]] = as.numeric(names(deg_possibilities[pos[1]:i]))
      }
      next

    }
    new_deg_list[[c]] = as.vector(unlist(deg_possibilities[pos[1]:i]))
    a = names(deg_possibilities[pos[1]]) # first degree
    b = names(deg_possibilities[(i)]) # last degree range
    names_deg[c] = paste(a,b,sep = "_to_")
    poss[[c]] = as.numeric(names(deg_possibilities[pos[1]:i]))

    c = c + 1

    times = 0
    l_sum = 0


    pos = numeric()
    p = 1

  }

  names(new_deg_list) = names_deg
  l_new = as.vector(sapply(new_deg_list, function(x) length(x)))

  # Mapping between degrees and where in the new_deg_list are located. So we know from where do I need to extract
  l_pos = as.vector(sapply(poss, function(x) length(x)))

  i = 1
  sub1 = numeric()
  sub = numeric()
  for(i in 1:length(l_pos)){

    sub = c(sub,poss[[i]])
    ok = rep(i,length(poss[[i]]))
    sub1 = c(sub1,ok)
  }

  # mapping done

  map_pos_deg_list = as.data.frame(cbind(sub, sub1))


  # geneset = data[[1]]

  # Take the query gene sets and generate randomizations from it. Input gene set and number of randomizations
  #Opposite to %in%
  '%!in%' <- function(x,y)!('%in%'(x,y))
  sampling_generator = function(geneset, n_randomizations) {
    b = as.data.frame(table(as.character(as.vector(geneset[[2]]))))
    b[, 1] = as.numeric(as.vector(b[, 1]))
    b = b[with(b, order(Var1)), ]
    gens = as.vector(geneset[[1]])
    a = lapply(new_deg_list, function(x) x[which(x %!in% gens)])
    # a = new_deg_list
    pos = map_pos_deg_list[which(map_pos_deg_list[, 1] %in%
                                   as.character(b[, 1])), 2]
    times = c(1:n_randomizations)


    ok = as.data.frame(cbind(pos,b[,2]))
    pp = unique(as.vector(ok[,1]))

    newlist = list()
    i = 1
    for(i in 1:length(pp)){

      ss1 =  sum(ok[which(ok[,1]%in%pp[i]),2])
      ss2 = pp[i]

      newlist[[i]] = c(ss2,ss1)
    }
    sub = do.call(rbind.data.frame, newlist)




    r_genes = list()
    i = 1
    for(i in 1:sampling){



      r_genes[[i]] = as.vector(unlist(mapply(FUN = function(x,
                                                            y) {

        set.seed(times[i])
        # set.seed(sample(times,1,replace = F))
        sample(x, y, replace = F)
      }, a[sub[,1]], sub[, 2], SIMPLIFY = FALSE)))

    }


    return(r_genes)
  }


  links_geneset = function(geneset) {
    links = links_matrix[geneset,]
    links = as.data.frame(links)
    links = sapply(links, function(x)
      sum(x))
    links = unname(links, force = FALSE)
    return(links)
  }



  index_genes = rownames(links_matrix)

  filtered_genesets = genesets[which(as.vector(genesets[,
                                                        1]) %in% as.vector(degreeList[,1])),]

  group_sets = unique(as.vector(filtered_genesets[, 2]))

  length_genesets_filtered = numeric()
  i = 1
  for (i in 1:length(group_sets)) {
    sub = filtered_genesets[which(filtered_genesets[, 2] %in% group_sets[i]),]
    length_genesets_filtered[i] = nrow(sub)
  }
  if (length(length_genesets_filtered) > 20) {
    set = seq(from = 1,
              to = length(length_genesets_filtered),
              by = 20)
    set[length(set)] = length(length_genesets_filtered)
  }else {
    set = length(length_genesets_filtered)
  }

  i = 1
  deg = numeric()
  for(i in 1:nrow(filtered_genesets)){

    deg[i] = as.numeric(as.vector(degreeList[which(degreeList[,
                                                              1] %in% as.vector(filtered_genesets[i, 1])), 2]))
  }

  filtered_genesets[,3] = deg

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
    for(i in 1:length(group_sets)){

      geneset_test = filtered_genesets[which(filtered_genesets[,2]%in%group_sets[i]),c(1,3)]
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


    # For the overlap
    i = 1
    overlap_list = list()
    for (i in 1:length(groups)) {
      sub = as.vector(genesets[which(genesets[,
                                              2] %in% groups[i]), 1])

      overlap_list[[i]] = sub
    }

    links_geneset_real = function(geneset) {
      links = links_matrix[geneset,]
      links = as.data.frame(links)
      links = sapply(links, function(x)
        sum(x))
      links = unname(links, force = FALSE)
      return(links)
    }
    links_geneset = function (geneset){


      links = links_matrix[which(rownames(links_matrix)%in%as.vector(unlist(geneset))),]
      links = as.data.frame(links)
      links = sapply(links, function(x) sum(x))
      links = unname(links, force = FALSE)
      return(links)
    }


    if (website == TRUE) {

      loglik = function(inits) {
        A <- inits[1]
        B <- inits[2]
        Y = dat[, 3]
        N = dat[, 2]
        - sum(
          lgamma(abs(A) + abs(B)) - lgamma(abs(A)) - lgamma(abs(B)) + lgamma(Y + abs(A)) + lgamma(N -
                                                                                                    Y + abs(B)) - lgamma(N + abs(A) + abs(B))
        )
      }

      total = 0
      result_t = data.frame()
      timer = 1
      cc = 1

      for (timer in 1:length(set)) {
        no_cores <- cores
        cl <- makeCluster(no_cores)


        times = sampling
        geneset_test_list1 = geneset_test_list[cc:set[timer]]
        data = parLapply(cl,geneset_test_list1, function(x) as.data.frame(x))
        clusterExport(cl,
                      list("sampling_generator", "new_deg_list","%!in%","map_pos_deg_list","times"),
                      envir = environment())
        length_geneset1 = length_genesets[cc:set[timer]]
        group_sets1 = group_sets[cc:set[timer]]

        prueba1 = parSapply(cl,data, function(x,y) sampling_generator(x,y), y = times)
        prueba1 = lapply(prueba1, function(x){as.vector(unlist(x))})

        clusterExport(cl, list("links_geneset", "links_matrix"), envir = environment())
        query = parLapply(cl, prueba1, function(x)
          links_geneset(x))

        m <- length(query[[1]])
        clusterExport(cl, list("m", "query"), envir = environment())
        remove(prueba1)
        information_list_true = parLapply(cl, 1:m, function(j)
          sapply(query, "[[", j))
        remove(query)
        stopCluster(cl)

        filtered_genesets = genesets[which(as.vector(genesets[, 1]) %in% index_genes),]
        groups = as.vector(unique(filtered_genesets[, 2]))
        i = 1
        query_list = list()
        for (i in 1:length(groups)) {
          sub = filtered_genesets[which(filtered_genesets[, 2] %in% groups[i]), 1]
          pos = which(index_genes %in% sub)
          query_list[[i]] = pos

        }
        real_genesets = lapply(query_list[cc:set[timer]], function(x)
          links_geneset_real(x))




        nose = numeric()
        names(real_genesets) = groups[cc:set[timer]]
        paths = as.vector(unique(pathways[, 2]))

        pvalue = numeric()
        exp_mean = numeric()
        gene_sets = vector()
        pathwys = vector()
        obv_links = numeric()
        overlapp = numeric()
        i = 1
        y = 1
        k = 1
        d = times
        for (i in 1:length(real_genesets)) {
          j = 1
          for (j in 1:length(information_list_true)) {
            obv = as.numeric(as.vector(real_genesets[[i]][j]))
            subset = information_list_true[[j]][y:d]
            g_set = as.vector(genesets[, 1][which(genesets[, 2] %in% names(real_genesets[i]))])
            overlap = length(g_set[which(g_set %in% as.vector(pathways[which(pathways[, 2] %in%
                                                                               paths[j]), 1]))])
            # if(is.na(overlap)==TRUE){overlap = 0}
            max = as.numeric((length_geneset1[i] * length_pathways[j]) - min(length_geneset1[i],length_pathways[j]))


            # if (mean(subset) == 0) {
            #   subset[times] = 1
            # }


            n = max
            m_1 = mean(subset)
            m_2 = mean(subset ^ 2)

            alpha = (n * m_1 - m_2) / (n * (m_2 / m_1 - m_1 - 1) + m_1)
            beta = (n - m_1) * (n - m_2 / m_1) / (n * (m_2 / m_1 - m_1 - 1) +
                                                    m_1)

            inits = c(alpha, beta)


            if (any(is.nan(inits)) == TRUE) {
              pvalue[k]= NA

            }else if (any(inits == 0) == TRUE) {
              pvalue[k] = NA

            }else{
              dat = as.data.frame(cbind(1:times, rep(max, times), subset))
              optim.tas = optimr(inits, loglik)

              optim.tas$par = abs(optim.tas$par)



              pvalue[k] = 0.5 * dbb(obv, max, optim.tas$par[1], optim.tas$par[2]) + sum(dbb((obv +
                                                                                               1):(max),
                                                                                            (max),
                                                                                            optim.tas$par[1],
                                                                                            optim.tas$par[2]))

            }



            exp_mean[k] = mean(subset)
            obv_links[k] = obv
            gene_sets[k] = names(real_genesets[i])
            pathwys[k] = paths[j]
            overlapp[k] = overlap
            k = k + 1

            # if callback is set, use it to send back the message:
            if (!is.null(callback)) {
              callback(paste(length(paths) - j, " pathways remaining", sep = ""))
            }
            print(paste(length(paths) - j, " pathways remaining", sep = ""))
          }
          y = y + times
          d = d + times
          total = total + 1
          if (!is.null(callback)) {
            callback(paste(sum(set) - total, " gene sets remaining", sep = ""))
          }
          print(paste(sum(set) - total, " gene sets remaining", sep = ""))
        }

        result = as.data.frame(cbind(gene_sets, pathwys, obv_links, exp_mean, overlapp, pvalue))
        result[, 7] = p.adjust(pvalue, method = "BH")
        # result[, 7] = p.adjust(pvalue, method = "BH")
        # result[, 7] = p.adjust(result[, 6], method = "bonferroni")
        result[, 8] = p.adjust(pvalue, method = "bonferroni")
        result[, 3] = as.numeric(as.vector(result[, 3]))
        result[, 4] = as.numeric(as.vector(result[, 4]))
        result[, 5] = as.numeric(as.vector(result[, 5]))
        result[, 6] = as.numeric(as.vector(result[, 6]))

        result_t = rbind(result_t, result)

        cc = set[timer] + 1



      }
      names(result_t) = c("geneset",
                          "pathway",
                          "obv_links",
                          "exp_mean",
                          "overlap",
                          "p-value",
                          "q-value",
                          "FWER")

    } else{
      total = 0
      timer = 1
      cc = 1
      result_t <- list()

      loglik = function(inits, x) {
        A <- inits[1]
        B <- inits[2]
        Y = x[, 3]
        N = x[, 2]
        - sum(
          lgamma(abs(A) + abs(B)) - lgamma(abs(A)) - lgamma(abs(B)) +
            lgamma(Y + abs(A)) + lgamma(N - Y + abs(B)) - lgamma(N +
                                                                   abs(A) + abs(B))
        )
      }


      stat_computation = function(c, observed1, max1) {
        n = max1
        m_1 = mean(sol[[c]])
        m_2 = mean(sol[[c]] ^ 2)
        alpha = (n * m_1 - m_2) / (n * (m_2 / m_1 - m_1 - 1) + m_1)
        beta = (n - m_1) * (n - m_2 / m_1) / (n * (m_2 / m_1 - m_1 -
                                                     1) + m_1)
        inits = c(alpha, beta)
        if (any(is.nan(inits)) == TRUE) {
          pvalue = NA
        }
        else if (any(inits == 0) == TRUE) {
          pvalue = NA
        }
        else {
          dat = as.data.frame(cbind(1:times, rep(max1, times),
                                    sol[[c]]))
          optim.tas = optimr(inits, fn = loglik, x = dat)
          optim.tas$par = abs(optim.tas$par)
          pvalue = 0.5 * dbb(observed1, max1, optim.tas$par[1],
                             optim.tas$par[2]) + sum(dbb((observed1 + 1):(max1),
                                                         (max1),
                                                         optim.tas$par[1],
                                                         optim.tas$par[2]
                             ))
        }
        return(pvalue)
      }


      for (timer in 1:length(set)) {

        no_cores <- cores
        cl <- makeCluster(no_cores)
        times = sampling
        geneset_test_list1 = geneset_test_list[cc:set[timer]]
        data = parLapply(cl,geneset_test_list1, function(x) as.data.frame(x))
        clusterExport(cl,
                      list("sampling_generator", "%!in%","new_deg_list","map_pos_deg_list","times"),
                      envir = environment())
        length_geneset1 = length_genesets_filtered[cc:set[timer]]
        group_sets1 = group_sets[cc:set[timer]]

        prueba1 = parSapply(cl,data, function(x,y) sampling_generator(x,y), y = times)
        prueba1 = lapply(prueba1, function(x){as.vector(unlist(x))})
        clusterExport(cl, list("links_geneset", "links_matrix"),
                      envir = environment())
        query = parLapply(cl, prueba1, function(x)
          links_geneset(x))
        m <- length(query[[1]])
        clusterExport(cl, list("m", "query"), envir = environment())
        remove(prueba1)
        information_list_true = parLapply(cl, 1:m, function(j)
          sapply(query,
                 "[[", j))
        remove(query)
        #
        t1 = Sys.time()
        real_genesets = lapply(query_list[cc:set[timer]], function(x)
          links_geneset_real(x))
        stopCluster(cl)


        if (!is.null(callback)) {
          callback("50% of the process done")
        }




        no_cores <- cores
        cl <- makeCluster(no_cores)
        clusterExport(
          cl,
          list(
            "information_list_true",
            "real_genesets",
            "query_list",
            "links_matrix",
            "times",
            "group_paths",
            "pathways",
            "length_pathways",
            "length_geneset1"
          ),
          envir = environment()
        )

        sol = parSapply(cl, information_list_true, function(x)
          split(x,
                ceiling(seq_along(x) /
                          times)))
        clusterExport(cl, list("sol"), envir = environment())

        remove(information_list_true)
        observed = as.vector(unlist(real_genesets))
        g_set = parLapply(cl, overlap_list[cc:set[timer]], function(x)
          as.vector(unlist(x)))

        overlapp = parLapply(cl, g_set, function(x)
          as.vector(sapply(group_paths,
                           function(y)
                             length(x[which(x %in% as.vector(pathways[which(pathways[,
                                                                                     2] %in% y), 1]))]))))
        overlapp = unlist(overlapp)
        remove(g_set)
        length_pathways = as.numeric(as.vector(length_pathways))
        length_genesets1 = as.numeric(as.vector(length_geneset1))
        max = unlist(parLapply(cl, length_geneset1, function(x)
          as.vector(
            sapply(length_pathways,
                   function(y)
                     (y * x) - min(y,x))
          )))
        # max = max - overlapp

        stopCluster(cl)


        if (!is.null(callback)) {
          callback("65% of the process done")
        }

        positions = seq(
          from = 1,
          to = length(length_pathways) * length(length_geneset1),
          by = length(length_geneset1)
        )
        positions = rep(positions, length(length_geneset1))
        summ = rep(0:(length(length_geneset1) - 1), each = length(length_pathways))
        positions = positions + summ
        n.cores = cores
        if (.Platform$OS.type == "windows") {
          n.cores = 1
        }

        pvalues = as.vector(unlist(
          mcmapply(
            function(x, y, z)
              stat_computation(x,
                               y, z),
            positions,
            observed,
            max,
            mc.cores = n.cores
          )
        ))
        if (!is.null(callback)) {
          callback("90% of the process done")
        }
        expected = as.numeric(unlist(lapply(sol, function(x)
          mean(x))))
        remove(sol)
        expected = expected[positions]
        result = as.data.frame(cbind(
          rep(group_sets1, each = length(length_pathways)),
          rep(group_paths, length(length_genesets1)),
          observed,
          expected,
          overlapp,
          pvalues
        ))

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
    }

    return(result_t)

  }



}
