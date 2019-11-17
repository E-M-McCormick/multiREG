#' Find individual-level solutions using final patterns. 
#' @param subdata Data from previous step.
#' @param yvarnames Endogenous variable names.
#' @param interactnames Names of interaction variables.
#' @param grppen Group path penalties.
#' @param output Output object
#' @param verbose Logical. If TRUE, algorithm will print progress to console.
#' @return Returns final paths for all individuals. 
#' @keywords internal  
ind_search = function(subdata,
                      yvarnames,
                      interactnames,
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
  
  #### Loop through Subjects Data for Individual Search ####
  for (sub in names(subdata)){
    if(verbose){print(paste0('Building individual-level model for ', sub, '.'), quote = FALSE)}
    tempdata = subdata[[sub]]
    for (varname in yvarnames){
      subset_index = !(colnames(tempdata) %in% varname |
                         colnames(tempdata) %in% interactnames[grepl(paste0('\\<',varname,'_by'), interactnames)] | 
                         colnames(tempdata) %in% interactnames[grepl(paste0('by_',varname,'\\>'), interactnames)] |
                         colnames(tempdata) %in% interactnames[grepl(paste0('by_',varname,'_by'), interactnames)])
      
      final_coefs = model_selection(x = as.matrix(tempdata[, subset_index]),
                                    y = tempdata[, colnames(tempdata) %in% varname],
                                    selection_crit = model_crit,
                                    alpha = alpha,
                                    penalty.factor = initial_penalties[subset_index, varname])
      
      for (predictor in rownames(final_coefs)[!rownames(final_coefs) %in% '(Intercept)']){
        if (final_coefs[predictor,] != 0){
          finalpaths[predictor, varname, sub] = final_coefs[predictor,]
        }
      }
    }
  }
  return(finalpaths)
}
