#' multiREG Check Variability
#' @return Returns variability check.
#' @param data raw data files
#' @keywords internal 
check_variability = function(data = NULL){
  variability = list()
  variability[['flag']] = FALSE
  for (sub in names(data)){
    temp = apply(data[[sub]], 2, function(x) { range(x, na.rm = TRUE) })
    isEqual = apply(temp, 2, function(x) { isTRUE(all.equal(x[1], x[2], tol = .005)) } )
    if (any(isEqual == TRUE)){
      variability[['flag']] = TRUE
    }
    variability[[sub]] = isEqual
  }
  return(variability)
}