#' Find subgroup solutions using final patterns. 
#' @return Returns subgroup membership, similarity matrix, modularity, and number of subgroups.
#' @keywords internal  
subgroup_search <- function(subdata,
                            indpaths,
                            sub_method){

print(paste0('Search for Subgroups'), quote = FALSE)

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

g            <- graph.adjacency(sim, mode = "undirected", weighted = TRUE)
weights      <- E(g)$weight
  
if (sub_method == "Walktrap")
    res        <- cluster_walktrap(g, weights = weights, steps = 4)
  
if (sub_method == "Infomap")
    res        <- cluster_infomap(g, e.weights = weights)
  
if (sub_method == "Edge Betweenness")
    res        <- cluster_edge_betweenness(g, weights = weights)
  
if (sub_method == "Fast Greedy")
    res        <- cluster_fast_greedy(g, weights = weights)
  
if (sub_method == "Label Prop")
    res        <- cluster_label_prop(g, weights = weights)
  
if (sub_method == "Leading Eigen")
    res        <- cluster_leading_eigen(g, weights = weights)
  
if (sub_method == "Louvain")
    res        <- cluster_louvain(g, weights = weights)
  
if (sub_method == "Spinglass")
    res        <- cluster_spinglass(g, weights = weights)

subgroup_results     <- list()
subgroup_results$sub_mem    <- data.frame(names      = names(membership(res)), 
                           sub_membership = as.numeric(membership(res)))
subgroup_results$sim         <- sim
subgroup_results$n_subgroups <- length(unique(na.omit(subgroup_results$sub_mem$sub_membership))) 
subgroup_results$modularity  <- modularity(res)
  
return(subgroup_results)
} 
  


