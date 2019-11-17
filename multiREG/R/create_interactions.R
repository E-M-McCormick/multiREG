#' Create Interactions
#' @param endog Matrix of endogenous variables (at time and lagged).
#' @param exog Matrix of exogenous variables (at time and lagged).
#' @param interactions List of specified interactions.
#' @return Returns subject data with new interactions and centered main effects if desired.
#' @keywords internal 
create_interactions = function(endog = NULL, exog = NULL, interactions = NULL){
  #### If no interactions specified, will skip everything else. ####
  if (is.null(interactions)){
    new_vars = cbind(endog, exog)
    return(new_vars)
  }
  
  #### If 'all' is present, will skip everything else. ####
  if (any(interactions == 'all')){
    tempdat = scale(cbind(endog, exog), center = TRUE, scale = FALSE)
    int_vars = vector()
    newnames = vector()
    for (i in 1:(ncol(tempdat))){for(j in i:ncol(tempdat)){
      int_vars = cbind(int_vars, tempdat[,i] * tempdat[,j])
      newnames = cbind(newnames, paste0(colnames(tempdat)[i], '_by_', colnames(tempdat)[j]))
    }}
    colnames(int_vars) = newnames
    new_vars = cbind(tempdat, int_vars)
    return(new_vars)
  }
  
  #### If 'all_cross' is present, will skip everything else. ####
  if (any(interactions == 'all_cross')){
    tempdat = scale(cbind(endog, exog), center = TRUE, scale = FALSE)
    int_vars = vector()
    newnames = vector()
    for (i in 1:(ncol(tempdat)-1)){for(j in (i+1):ncol(tempdat)){
      int_vars = cbind(int_vars, tempdat[,i] * tempdat[,j])
      newnames = cbind(newnames, paste0(colnames(tempdat)[i], '_by_', colnames(tempdat)[j]))
    }}
    colnames(int_vars) = newnames
    new_vars = cbind(tempdat, int_vars)
    return(new_vars)
  }
  
  #### Loop through desired interactions and create while centering any variable that is used to create an interaction. ####
  interact_vars = vector()
  for (q in 1:length(interactions)){
    int_vars = vector()
    newnames = vector()
    if (interactions[[q]] == 'all_exogenous'){
      exog = scale(exog, center = TRUE, scale = FALSE)
      for (i in 1:(ncol(exog)-1)){for(j in (i+1):ncol(exog)){
        int_vars = cbind(int_vars, exog[,i] * exog[,j])
        newnames = cbind(newnames, paste0(colnames(exog)[i], '_by_', colnames(exog)[j]))
      }}
      colnames(int_vars) = newnames
    } else if (interactions[[q]] == 'all_endogenous'){
      endog = scale(endog, center = TRUE, scale = FALSE)
      for (i in 1:(ncol(endog)-1)){for(j in (i+1):ncol(endog)){
        int_vars = cbind(int_vars, endog[,i] * endog[,j])
        newnames = cbind(newnames, paste0(colnames(endog)[i], '_by_', colnames(endog)[j]))
      }}
      colnames(int_vars) = newnames
    } else if (interactions[[q]] == 'all_endog_by_exog'){
      endog = scale(endog, center = TRUE, scale = FALSE)
      exog = scale(exog, center = TRUE, scale = FALSE)
      for (i in 1:(ncol(endog))){for(j in 1:ncol(exog)){
        int_vars = cbind(int_vars, endog[,i] * exog[,j])
        newnames = cbind(newnames, paste0(colnames(endog)[i], '_by_', colnames(exog)[j]))
      }}
      colnames(int_vars) = newnames
    } else {
      int_splitnames = unlist(strsplit(interactions[[q]] ,'\\*'))
      endog[, colnames(endog) %in% int_splitnames] = scale(endog[, colnames(endog) %in% int_splitnames, drop = FALSE], center = TRUE, scale = FALSE)
      exog[, colnames(exog) %in% int_splitnames] = scale(exog[, colnames(exog) %in% int_splitnames, drop = FALSE], center = TRUE, scale = FALSE)
      tempmerg = vector()
      for (r in 1:length(int_splitnames)){
        tempmerg = cbind(tempmerg, cbind(endog[, colnames(endog) %in% int_splitnames[r], drop = FALSE], exog[, colnames(exog) %in% int_splitnames[r], drop = FALSE]))
      }
      int_vars = as.matrix(apply(tempmerg, 1, prod))
      colnames(int_vars) = paste(int_splitnames, collapse = '_by_')
    }
    interact_vars = cbind(interact_vars, int_vars)
  }
  
  #### Check for and Remove Duplicated Interaction Vars ####
  interact_vars = interact_vars[, !duplicated(round(interact_vars, 10), MARGIN = 2), drop = FALSE]
  
  new_vars = cbind(endog, exog, interact_vars)
  return(new_vars)
}
