#' Find subgroup solutions using final patterns. 
#' @return Returns subgroup membership, similarity matrix, modularity, and number of subgroups.
#' @keywords internal  
group_search <- function(subdata,
                         groupcutoff,
                         yvarnames,
                         interact_exogenous,
                         predict_with_interactions,
                         interactnames,
                         interact_exogvars, 
                         output,
                         grppen = NULL,
                         initial_penalties = NULL){
  
  numvars = ncol(subdata[[1]])
  model_crit = output[['function_parameters']][['model_crit']]
  alpha = output[['function_parameters']][['alpha']]

  pathpresent = group_coefs = array(data = rep(NaN, numvars*numvars*length(subdata)), 
                                    dim = c(numvars, numvars, length(subdata)), 
                                    dimnames = list(c(colnames(subdata[[1]])), 
                                                    c(colnames(subdata[[1]])), 
                                                    c(names(subdata))))
# Loop through Subjects Data for Group Search
for (sub in names(subdata)){
  if(is.null(grppen))
    print(paste0('Building group-level model for ', sub, '.'), quote = FALSE)
  if(!is.null(grppen))
    print(paste0('Building subgroup-level model for ', sub, '.'), quote = FALSE)
  tempdata = subdata[[sub]]
  for (varname in yvarnames){
    subset_predictors = as.matrix(tempdata[, !(colnames(tempdata) %in% varname |
                                                 colnames(tempdata) %in% paste0(varname,'_by_',interact_exogvars))])
    if (!is.null(predict_with_interactions) & !varname %in% predict_with_interactions){
      subset_predictors = subset_predictors[, !colnames(subset_predictors) %in% interactnames]
    }
    
    if(is.null(grppen))
    final_coefs = model_selection(x = subset_predictors,
                                  y = tempdata[, colnames(tempdata) %in% varname],
                                  selection_crit = model_crit,
                                  alpha = alpha,
                                  penalty.factor = initial_penalties[!colnames(tempdata) %in% varname, varname])
    # if subgroup stage, use grppen to start
    #set up things if subgroup-level stage (indicated via use of grppen)
    if(!is.null(grppen)){
      subset_penalties = grppen[!(colnames(tempdata) %in% varname |
                                    colnames(tempdata) %in% paste0(varname,'_by_',interact_exogvars)), varname] 
      if (!is.null(predict_with_interactions) & !varname %in% predict_with_interactions){
        subset_penalties = subset_penalties[!names(subset_penalties) %in% interactnames] 
      }
    if(!is.null(grppen))
    final_coefs = model_selection(x = subset_predictors,
                                  y = tempdata[, colnames(tempdata) %in% varname],
                                  selection_crit = model_crit,
                                  alpha = alpha,
                                  penalty.factor = subset_penalties)
    }
    for (predictor in rownames(final_coefs)[!rownames(final_coefs) %in% '(Intercept)']){
      if (final_coefs[predictor,] == 0){
        pathpresent[predictor, varname, sub] = 0
      } else {
        pathpresent[predictor, varname, sub] = 1
      }
      group_coefs[predictor, varname, sub] = final_coefs[predictor, ]
    }
  }
}

# Calculate Paths that Should Appear in the Group (Non-Penalized) Model
group_thresh_mat = rowSums(pathpresent, dims = 2)
group_thresh_mat = group_thresh_mat/(length(subdata))
group_thresh_mat[group_thresh_mat < groupcutoff] = 0
group_thresh_mat[group_thresh_mat >= groupcutoff] = 1
group_penalties = abs(group_thresh_mat - 1)
output[['group']][['group_paths_counts']] = group_thresh_mat
output[['group']][['group_paths_proportions']] = group_thresh_mat

grppaths <- list()
grppaths <- list("output" = output, 
                 "group_thresh_mat" = group_thresh_mat, 
                 "group_penalties" = group_penalties)

return(grppaths)
}
