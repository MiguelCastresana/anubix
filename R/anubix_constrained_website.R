
#' Perform ANUBIX for the provided gene sets
#'
#' @description Computes ANUBIX, an accurate test for network enrichment analysis between query sets and pathway sets. Instead of normal random sampling it does constrained random sampling, taking the degree of the nodes into account.
#' @param links_matrix A numeric matrix. Links_matrix stores the links each gene has to each pathway. The number of rows equals to the number of genes in the network and the number of columns equals to the number of pathways under study.
#' @param genesets A data.frame of two columns. Column 1, are the genes and column 2 the experiment where they belong to.
#' @param pathways A data.frame of two columns. Column 1 are the genes and second column the pathway where they belong to.
#' @param new_deg_list A list with the genes that fall into certain degree
#' @param map_pos_deg_list A data frame which maps the position of a degree value in the new_deg_list list
#' @param cores A numeric value. Cores is defined as the number of cores used by the algorithm. The default value is \strong{2}
#' @param sampling A numeric value. Sampling is defined as the number of random samplings required to construct the null distribution. The default value is \strong{2000}
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
#' @seealso \code{\link{anubix_links}},\code{\link{example_anubix}},\code{\link{anubix_clustering}},\code{\link{anubix}}
#'





anubix_constrained_website = function(links_matrix,
                              genesets,
                              pathways,
                              cores = 2, 
                              sampling = 2000,
                              callback = NULL) {
  

    

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
  
  
  
  # Take the query gene sets and generate randomizations from it. Input gene set and number of randomizations
  #Opposite to %in%
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  
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
    for(i in 1:2000){
      
      
      
      r_genes[[i]] = as.vector(unlist(mapply(FUN = function(x, 
                                                            y) {
        
        set.seed(times[i])
        # set.seed(sample(times,1,replace = F))
        sample(x, y, replace = F)
      }, a[sub[,1]], sub[, 2], SIMPLIFY = FALSE)))
      
    }
    
    
    return(r_genes)
  }
  
  

  set = length(length_genesets)
  

  
  
  
  index_genes = rownames(links_matrix)
  
  filtered_genesets = genesets[which(as.vector(genesets[,
                                                        1]) %in% as.vector(degreeList[,1])),]
  
  
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
    
    group_sets = unique(as.vector(filtered_genesets[, 2]))
    
    length_genesets_filtered = numeric()
    i = 1
    for (i in 1:length(group_sets)) {
      sub = filtered_genesets[which(filtered_genesets[, 2] %in% group_sets[i]),]
      length_genesets_filtered[i] = nrow(sub)
    }
    
    
    i = 1
    deg = numeric()
    for(i in 1:nrow(filtered_genesets)){
      
      deg[i] = as.numeric(as.vector(degreeList[which(degreeList[, 
                                                                1] %in% as.vector(filtered_genesets[i, 1])), 2]))
    }
    
    filtered_genesets[,3] = deg
    
    
    
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
      links = sapply(links, function(x)
        sum(x))
      links = unname(links, force = FALSE)
      return(links)
    }
    links_geneset = function (geneset){
      
      
      links = links_matrix[which(rownames(links_matrix)%in%as.vector(unlist(geneset))),]
      links = sapply(links, function(x) sum(x))
      links = unname(links, force = FALSE)
      return(links)
    }
    
    
    loglik = function(inits) {
      A <- inits[1]
      B <- inits[2]
      - sum(
        lgamma(abs(A) + abs(B))
        - lgamma(abs(A))
        - lgamma(abs(B))
        + lgamma(Y + abs(A))
        + lgamma(N - Y + abs(B))
        - lgamma(n0 + abs(A) + abs(B))
      )
    }
    
    
    
    
    no_cores <- cores
    cl <- makeCluster(no_cores)
    
    
    times = sampling
    geneset_test_list1 = geneset_test_list
    data = parLapply(cl,geneset_test_list1, function(x) as.data.frame(x))
    clusterExport(cl,
                  list("sampling_generator", "new_deg_list","%!in%","map_pos_deg_list","times"),
                  envir = environment())
    length_geneset1 = length_genesets
    group_sets1 = group_sets
    
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
    real_genesets = lapply(query_list, function(x)
      links_geneset_real(x))
    
    
    
    total = 0
    
    names(real_genesets) = groups
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
        
        max = as.numeric((length_geneset1[i] * length_pathways[j]) - min(length_geneset1[i],length_pathways[j]))
        
        
        
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
          Y = dat[, 3]
          N = dat[, 2]
          n0 = dat[1, 2]
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
    
    result[, 8] = p.adjust(pvalue, method = "bonferroni")
    result[, 3] = as.numeric(as.vector(result[, 3]))
    result[, 4] = as.numeric(as.vector(result[, 4]))
    result[, 5] = as.numeric(as.vector(result[, 5]))
    result[, 6] = as.numeric(as.vector(result[, 6]))
    
    
    
    names(result) = c("geneset",
                      "pathway",
                      "obv_links",
                      "exp_mean",
                      "overlap",
                      "p-value",
                      "q-value",
                      "FWER")
    
    
   
  
    
  }
  
  return(result)

}
