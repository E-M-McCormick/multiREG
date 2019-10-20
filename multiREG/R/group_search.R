#' Find subgroup solutions using final patterns.
#' @param subdata Data from previous step.
#' @param groupcutoff Thresholds for group-level paths.
#' @param yvarnames Endogenous variable names.
#' @param interact_exogenous Exogenous variable to use as moderators.
#' @param predict_with_interactions Endogenous variables to be predicted by interactions.
#' @param interactnames Names of interaction variables.
#' @param interact_exogvars Exogenous variable to use as moderators.
#' @param output Output object
#' @param grppen Group path penalties.
#' @param initial_penalties User-set penalties.
#' @param verbose Logical. If TRUE, algorithm will print progress to console.
#' @return Returns subgroup membership, similarity matrix, modularity, and number of subgroups.
#' @keywords internal  
group_search = function(subdata,
                        groupcutoff,
                        yvarnames,
                        interact_exogenous,
                        predict_with_interactions,
                        interactnames,
                        interact_exogvars, 
                        output,
                        grppen = NULL,
                        initial_penalties = NULL,
                        verbose){
  
  numvars = ncol(subdata[[1]])
  model_crit = output[['function_parameters']][['model_crit']]
  alpha = output[['function_parameters']][['alpha']]

  pathpresent = array(data = rep(NaN, numvars*numvars*length(subdata)), 
                      dim = c(numvars, numvars, length(subdata)), 
                      dimnames = list(c(colnames(subdata[[1]])), 
                                      c(colnames(subdata[[1]])), 
                                      c(names(subdata))))
  
# Loop through Subjects Data for Group Search
for (sub in names(subdata)){
  if(is.null(grppen))
    if(verbose){print(paste0('Building group-level model for ', sub, '.'), quote = FALSE)}
  if(!is.null(grppen))
    if(verbose){print(paste0('Building subgroup-level model for ', sub, '.'), quote = FALSE)}
  tempdata = subdata[[sub]]
  for (varname in yvarnames){
    subset_predictors = as.matrix(tempdata[, !(colnames(tempdata) %in% varname |
                                                 colnames(tempdata) %in% paste0(varname,'_by_',interact_exogvars))])
    if (!is.null(predict_with_interactions) & !varname %in% predict_with_interactions){
      subset_predictors = subset_predictors[, !colnames(subset_predictors) %in% interactnames]
    }
    
    # For Initial Group Search
    if(is.null(grppen)){
      final_coefs = model_selection(x = subset_predictors,
                                    y = tempdata[, colnames(tempdata) %in% varname],
                                    selection_crit = model_crit,
                                    alpha = alpha,
                                    penalty.factor = initial_penalties[!colnames(tempdata) %in% varname, varname])
    }
    
    
    # For Subgroup Search (i.e., use grppen to start)
    if(!is.null(grppen)){
      subset_penalties = grppen[!(colnames(tempdata) %in% varname |
                                    colnames(tempdata) %in% paste0(varname,'_by_',interact_exogvars)), varname] 
      if (!is.null(predict_with_interactions) & !varname %in% predict_with_interactions){
        subset_penalties = subset_penalties[!names(subset_penalties) %in% interactnames] 
      }
      final_coefs = model_selection(x = subset_predictors,
                                    y = tempdata[, colnames(tempdata) %in% varname],
                                    selection_crit = model_crit,
                                    alpha = alpha,
                                    penalty.factor = subset_penalties)
    }
    
    # Indicate whether paths exist for this subject (for creating group penalty matrix)
    for (predictor in rownames(final_coefs)[!rownames(final_coefs) %in% '(Intercept)']){
      if (final_coefs[predictor,] == 0){
        pathpresent[predictor, varname, sub] = 0
      } else {
        pathpresent[predictor, varname, sub] = 1
      }
    }
  }
}

# Calculate Paths that Should Appear in the Group (Non-Penalized) Model
group_thresh_mat = rowSums(pathpresent, dims = 2)/(length(subdata))
group_thresh_mat[group_thresh_mat < groupcutoff] = 0
group_thresh_mat[group_thresh_mat >= groupcutoff] = 1
group_penalties = abs(group_thresh_mat - 1)

grppaths <- list()
grppaths <- list("output" = output, 
                 "group_thresh_mat" = group_thresh_mat, 
                 "group_penalties" = group_penalties)

return(grppaths)
}
