#' Find subgroup solutions using final patterns. 
#' @return Returns subgroup membership, similarity matrix, modularity, and number of subgroups.
#' @keywords internal  
subgroup_search <- function(subdata,
                            indpaths,
                            sub_method){

print(paste0('Searching for subgroups.'), quote = FALSE)

binary <- ifelse(abs(indpaths)> 0, 1, 0)

sim <- matrix(0, ncol = length(subdata), nrow = length(subdata))

for (sub in 1:length(subdata)){
  for (sub2 in 1:length(subdata)){
    sim[sub,sub2] <- sum(binary[,,sub] == 1 & binary[,,sub2] == 1 & 
                         sign(indpaths[,,sub]) == sign(indpaths[,,sub2]), na.rm = TRUE)
  }
}

sim <- sim-min(sim, na.rm = TRUE)

colnames(sim) <- rownames(sim) <- names(subdata)

g            <- igraph::graph.adjacency(sim, mode = "undirected", weighted = TRUE)
weights      <- igraph::E(g)$weight
  
if (sub_method == "Walktrap")
    res        <- igraph::cluster_walktrap(g, weights = weights, steps = 4)
  
if (sub_method == "Infomap")
    res        <- igraph::cluster_infomap(g, e.weights = weights)
  
if (sub_method == "Edge Betweenness")
    res        <- igraph::cluster_edge_betweenness(g, weights = weights)
  
if (sub_method == "Fast Greedy")
    res        <- igraph::cluster_fast_greedy(g, weights = weights)
  
if (sub_method == "Label Prop")
    res        <- igraph::cluster_label_prop(g, weights = weights)
  
if (sub_method == "Leading Eigen")
    res        <- igraph::cluster_leading_eigen(g, weights = weights)
  
if (sub_method == "Louvain")
    res        <- igraph::cluster_louvain(g, weights = weights)
  
if (sub_method == "Spinglass")
    res        <- igraph::cluster_spinglass(g, weights = weights)

subgroup_results     <- list()
subgroup_results$sub_mem    <- data.frame(names      = names(igraph::membership(res)), 
                           sub_membership = as.numeric(igraph::membership(res)))
subgroup_results$sim         <- sim
subgroup_results$n_subgroups <- length(unique(na.omit(subgroup_results$sub_mem$sub_membership))) 
subgroup_results$modularity  <- igraph::modularity(res)
  
return(subgroup_results)
} 
  


