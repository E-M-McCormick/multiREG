#' Find individual-level solutions using final patterns. 
#' @param subdata Data from previous step.
#' @param yvarnames Endogenous variable names.
#' @param interact_exogenous Exogenous variable to use as moderators.
#' @param predict_with_interactions Endogenous variables to be predicted by interactions.
#' @param interactnames Names of interaction variables.
#' @param interact_exogvars Exogenous variable to use as moderators.
#' @param grppen Group path penalties.
#' @param output Output object
#' @param verbose Logical. If TRUE, algorithm will print progress to console.
#' @return Returns final paths for all individuals. 
#' @keywords internal  
ind_search = function(subdata,
                      yvarnames,
                      interact_exogenous,
                      predict_with_interactions,
                      interactnames,
                      interact_exogvars,
                      grppen,
                      output,
                      verbose){
  
  numvars = ncol(subdata[[1]])
  model_crit = output[['function_parameters']][['model_crit']]
  alpha = output[['function_parameters']][['alpha']]
  
  finalpaths = array(data = rep(0, numvars*numvars*length(subdata)), 
                     dim = c(numvars, numvars, length(subdata)),
                     dimnames = list(c(colnames(subdata[[1]])),
                                     c(colnames(subdata[[1]])),
                                     c(names(subdata))))
  for (sub in names(subdata)){
    if(verbose){print(paste0('Building individual-level model for ', sub, '.'), quote = FALSE)}
    tempdata = subdata[[sub]]
    for (varname in yvarnames){
      subset_predictors = as.matrix(tempdata[, !(colnames(tempdata) %in% varname |
                                                   colnames(tempdata) %in% paste0(varname,'_by_',interact_exogvars))])
      if (!is.null(predict_with_interactions) & !varname %in% predict_with_interactions){
        subset_predictors = subset_predictors[, !colnames(subset_predictors) %in% interactnames]
      }
      
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
      
      for (predictor in rownames(final_coefs)[!rownames(final_coefs) %in% '(Intercept)']){
        if (final_coefs[predictor,] != 0){
          finalpaths[predictor, varname, sub] = final_coefs[predictor,]
        }
      }
    }
  }
  return(finalpaths)
}
