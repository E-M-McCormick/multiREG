#' Identify Moderated Group Paths
#' @param group_paths Matrix indicating which paths are group-level.
#' @param finalpaths Matrix of final coefficients for each path.
#' @param subgroup_mem Vector of subgroup membership.
#' @return Returns indices of moderated group paths.
#' @keywords internal
group_subgroupsign = function(group_paths,
                                 finalpaths,
                                 subgroup_mem){
  
  group_subgroupsign = group_paths
  group_subgroupsign[group_subgroupsign != 1] = 1
  groupindex = which(group_paths == 1, arr.ind=TRUE)
  
  for (index in 1:nrow(groupindex)){
    groupcoefs = finalpaths[groupindex[index,1], groupindex[index,2], ]
    groupmeans = stats::aggregate(groupcoefs, subgroup_mem, mean)
    if (any(sign(groupmeans$x) == 1 && sign(groupmeans$x) == -1)){
      amodel = stats::aov(groupcoefs ~ subgroup_mem)
      if (summary(amodel)[[1]][['Pr(>F)']][[1]] < .05){group_subgroupsign[groupindex[index,1], groupindex[index,2]] = -1}
    }
  }
  return(group_subgroupsign)
}