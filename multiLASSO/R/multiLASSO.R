#' @name multiLASSO
#' @aliases multiLASSO
#' @title Network Model Search using Regularization with LASSO
#' @description This function utilizes regression with regularization to build models for individuals 
#' consisting of individual and group-level paths.
#' @usage 
#' multiLASSO(data                       = NULL,
#'            out                        = NULL,
#'            sep                        = NULL,
#'            header                     = TRUE,
#'            ar                         = TRUE,
#'            plot                       = TRUE,
#'            conv_vars                  = NULL,
#'            conv_length                = 16,
#'            conv_interval              = 1,
#'            groupcutoff                = .75,
#'            alpha                      = .5,
#'            model_crit                 = 'bic',
#'            penalties                  = NULL,
#'            test_penalties             = FALSE,
#'            exogenous                  = NULL,
#'            lag_exogenous              = FALSE,
#'            interact_exogenous         = NULL,
#'            interact_with_exogenous    = NULL,
#'            predict_with_interactions  = NULL,
#'            subgroup                   = FALSE,
#'            sub_method                 = 'Walktrap')
#'             
#' @param data The path to the directory where individual data files are located,
#' or the name of the list containing individual data. Each file or matrix within the list
#' must contain a single matrix containing the a T (time) by p (number of variables) matrix,
#' where the rows represent time and columns represent individual variables. Individuals may
#' have different numbers of observations (T), but must have the same number of variables (p).
#' 
#' @param out (Optional) The path to directory where results will be stored. If specified,
#' a copy of output data will be saved into the directory. If the specified directory does
#' not exist, it will be created.
#' 
#' @param sep Spacing scheme for input files. 
#' '' indicates space-separated; ',' indicates comma separated; '/t' indicates tab-separated
#' Only necessary when reading in files from physical directory.
#' 
#' @param header (Logical) Indicate TRUE if variable names included in input file, FALSE otherwise.
#' Only necessary when reading in files from physical directory.
#' 
#' @param ar (Logical) If TRUE, begin model search with all autoregressive pathways estimated
#' with no shrinkage (i.e., penalty = 0).
#' 
#' @param plot (Logical) IF TRUE, will create pdf plots of network maps during output.
#' 
#' @param conv_vars Vector of variable names to be convolved via smoothed Finite Impulse 
#' Response (sFIR). Note, conv_vars are not not automatically considered exogenous variables.
#' To treat conv_vars as exogenous use the exogenous argument. Variables listed in conv_vars 
#' must be binary variables. If there is missing data in the endogenous variables their values 
#' will be imputed for the convolution operation only. Defaults to NULL. ### If there are multiple 
#' variables listed in conv_vars they are not used in the convolution of additional conv_vars.## 
#' You can't do lagged variables.
#' 
#' @param conv_length Expected response length in seconds. For functional MRI BOLD, 16 seconds (default) is typical
#' for the hemodynamic response function. 
#' 
#' @param conv_interval Interval between data acquisition. Currently must be a constant. For 
#' fMRI studies, this is the repetition time. Defaults to 1. 
#' 
#' @param groupcutoff Cutoff value for inclusion of a given path at the group-level.
#' For instance, group_cutoff = .75 indicates that a path needs to be estimated for 75% of
#' individuals to be included as a group-level path.
#' 
#' @param alpha Elastic-net parameter for the regularization approach. Values close to 0 mimic 
#' the ridge penalty, which tends to shrink correlated parameters towards one another. Values 
#' close to 1 mimic the lasso penalty, which tends to select one parameter and shrink
#' the others. The default value (alpha=.5) balances these two considerations, and tends to select
#' groups of correlated parameters and shrink other groups towards zero.
#' 
#' @param model_crit Argument to indicate the model selection criterion to use for model selection.
#' Defaults to 'bic' (select on BIC). Options: 'bic', 'aic', 'aicc', 'hqc', 'cv' (cross-validation).
#' 
#' @param penalties (Optional) A matrix of user-provided penalties to initialize group-model search. 
#' Should contain a column for all variables (including lagged versions and interactions) that will 
#' be included in the model search. Values of 1 (the default) will initialize a variable to be 
#' normally considered in the regularization, values of 0 will initialize a variable to be estimated
#' (i.e., no shrinkage), and values of Inf will exclude variables from the model.
#' 
#' @param test_penalties (Optional, Logical) Optional argument to output a sample penalty matrix
#' based on function parameters. Helpful for specifying a matrix to use in the penalties argument.
#' Function will exit gracefully before running anything if test_penalties = TRUE.
#' 
#' @param exogenous (Optional) A list of user-specified variables to consider as exogenous
#' (e.g., cannot be predicted) in the model search procedure. If variable names are supplied,
#' variables should be referred to by name. If not, then variables should be referenced by
#' the pattern 'V#', where # represents the column number in the original data file (e.g., 'V5').
#' 
#' @param lag_exogenous (Optional, Logical) If TRUE, a lagged version of the exogenous variable 
#' will be created. If set to TRUE, but exogenous variables are not indicated in the argument
#' above, the function will not run properly.
#'
#' @param interact_exogenous (Optional) Select which exogenous variables are used to create
#' interaction terms with other predictors (excluding other exogenous variables). If 
#' lag_exogenous = TRUE, then can specify both contemporaneous and lagged versions of exogenous
#' variables (e.g., c('V5', 'V4Lag')). Specify as 'all' to create interactions using all exogenous 
#' variables (contemporaneous or lagged) present in the data. Interactions between exogenous variables 
#' should be created manually and included as separate columns in the input file or list.
#' 
#' @param interact_with_exogenous (Optional) Select which endogenous variables are combined
#' with the exogenous variables specified with interact_exogenous to create interaction variables
#' (e.g., c('V1','V2','V2Lag')). Specify as 'all' to create interaction using all endogenous variables
#' in the data.
#' 
#' @param predict_with_interactions (Optional) Select which endogenous variables should be predicted
#' by interaction variables. This option cannot be used if interact_exogenous or interact_with_exogenous
#' are NULL.
#' 
#' @param subgroup Logical. If TRUE, subgroups are generated based on
#' similarities in model features using the \code{walktrap.community}
#' function from the \code{igraph} package. When ms_allow=TRUE, subgroup
#' should be set to FALSE.  Defaults to FALSE.  
#' 
#' @param sub_method Community detection method used to cluster individuals into subgroups. Options align 
#' with those available in the igraph package: "Walktrap" (default), "Infomap", "Louvain", "Edge Betweenness", 
#' "Label Prop", "Fast Greedy", "Leading Eigen", and "Spinglass". 
#' 
#' @import utils stats grDevices gimme
#' 
#' @export multiLASSO

multiLASSO = function(data                       = NULL,
                      out                        = NULL,
                      sep                        = NULL,
                      header                     = TRUE,
                      ar                         = TRUE,
                      plot                       = TRUE,
                      conv_vars                  = NULL,
                      conv_length                = 16,
                      conv_interval              = 1,
                      groupcutoff                = .75,
                      alpha                      = .5,
                      model_crit                 = 'bic',
                      penalties                  = NULL,
                      test_penalties             = FALSE,
                      exogenous                  = NULL,
                      lag_exogenous              = FALSE,
                      interact_exogenous         = NULL,
                      interact_with_exogenous    = NULL,
                      predict_with_interactions  = NULL,
                      subgroup                   = FALSE,
                      sub_method                 = 'Walktrap'){

  # Add Function Parameters to Output
  output = list()
  output[['function_parameters']] = as.list(sys.call())
  
  # Wrangle Data into List
  print('Reading in data.', quote = FALSE)
  if (!is.list(data)){
    subdata = list(); setwd(data)
    for (i in list.files(data)){
      tempname = tools::file_path_sans_ext(i)
      print(paste0('   Reading in ', tempname, '.'), quote = FALSE)
      subdata[[tempname]] = read.delim(i, sep=sep, header=header)
    } 
  } else if (is.list(data)){
      subdata = data
  }
  
  # Generate Variable Names as Needed
  varnames = colnames(subdata[[1]])
  if(is.null(varnames)){
    varnames = c(paste0('V', seq(1, length(subdata[[1]][1,]))))
    subdata = lapply(subdata, function(x){ colnames(x) = varnames; x })
  }

  # Convolve if indicated.
  if(!is.null(conv_vars)){
    varLabels <- list(
      conv = conv_vars, # variables to be convolved
      exog = exogenous, # user-specified exogenous variables
      coln = varnames   # all variable names
    )
    
    subdata <- gimme:::setupConvolve(
      ts_list       = subdata, 
      varLabels     = varLabels, 
      conv_length   = conv_length, 
      conv_interval = conv_interval
    )
  }
  
  # Categorize Variables. & Omit NaN Rows
  for (i in 1:length(subdata)){
    yvar = subdata[[i]][, !colnames(subdata[[i]]) %in% exogenous, drop=FALSE]
    yvarnames = colnames(yvar)
    lagvar = rbind(rep(NA, ncol(yvar)), yvar[1:(nrow(yvar)-1),])
    colnames(lagvar) = paste(colnames(lagvar), 'Lag', sep='')
    exogvar = subdata[[i]][, colnames(subdata[[i]]) %in% exogenous, drop=FALSE]
    
    # Set Endogenous to be Interacted
    if (!is.null(interact_with_exogenous) && interact_with_exogenous == 'all'){
      interact_with_exogenous = c(colnames(yvar), colnames(lagvar))
    }
    
    # Lag Exogenous Variables if Needed
    if (lag_exogenous == TRUE){
      lagexogvar = rbind(rep(NA, ncol(exogvar)), exogvar[1:(nrow(exogvar)-1), , drop=FALSE])
      colnames(lagexogvar) = paste(colnames(exogvar), 'Lag', sep='')
      if (!is.null(interact_exogenous) && interact_exogenous == 'all'){
        interact_exogvars = c(colnames(exogvar), colnames(lagexogvar))
      } else {
        interact_exogvars = interact_exogenous
      }
      # Create Interaction Variables As Needed
      if (!is.null(interact_exogenous)){
        
        tempendo = as.matrix(cbind(yvar[, colnames(yvar) %in% interact_with_exogenous, drop=FALSE], lagvar[,colnames(lagvar) %in% interact_with_exogenous, drop=FALSE]))
        tempexo = as.matrix(cbind(exogvar[, colnames(exogvar) %in% interact_exogvars, drop=FALSE], lagexogvar[,colnames(lagexogvar) %in% interact_exogvars, drop=FALSE]))
        
        interact_var = t(sapply(1:nrow(yvar), function(i){ tcrossprod( tempendo[i,], tempexo[i,]) }))
        
        newnames = vector()
        for (j in 1:ncol(tempexo)){ for (k in 1:ncol(tempendo)){
          newnames = append(newnames, paste0(colnames(tempendo)[k], '_by_', colnames(tempexo)[j]), length(newnames))
          } }
        
        colnames(interact_var) = interactnames = newnames
        subdata[[i]] = na.omit(cbind(yvar, lagvar, exogvar, lagexogvar, interact_var))
      } else {
        interact_exogvars = interactnames = ''
        subdata[[i]] = na.omit(cbind(yvar, lagvar, exogvar, lagexogvar))
      }
      exognames = c(colnames(exogvar), colnames(lagexogvar))
    } else if (lag_exogenous == FALSE){
      if (!is.null(interact_exogenous) && interact_exogenous == 'all'){
        interact_exogvars = colnames(exogvar)
      } else {
        interact_exogvars = interact_exogenous
      }
      # Create Interaction Variables As Needed
      if (!is.null(interact_exogenous)){
        
        tempendo = as.matrix(cbind(yvar[, colnames(yvar) %in% interact_with_exogenous, drop=FALSE], lagvar[, colnames(lagvar) %in% interact_with_exogenous, drop=FALSE]))
        tempexo = as.matrix(exogvar[, colnames(exogvar) %in% interact_exogvars, drop=FALSE])
        
        interact_var = t(sapply(1:nrow(yvar), function(i){ tcrossprod( tempendo[i,], tempexo[i,]) }))
        
        newnames = vector()
        for (j in 1:ncol(tempexo)){ for (k in 1:ncol(tempendo)){
          newnames = append(newnames, paste0(colnames(tempendo)[k], '_by_', colnames(tempexo)[j]), length(newnames))
        } }
        
        colnames(interact_var) = interactnames = newnames
        subdata[[i]] = na.omit(cbind(yvar, lagvar, exogvar, interact_var))
      } else {
        interact_exogvars = interactnames = ''
        subdata[[i]] = na.omit(cbind(yvar, lagvar, exogvar))
      }
      exognames = colnames(exogvar)
    }
  }
  print('Data successfully read in.', quote = FALSE)
  output[['variablenames']][['y_vars']] = yvarnames
  output[['variablenames']][['exogenous_vars']] = exognames
  output[['variablenames']][['interaction_vars']] = interactnames
  
  # Check for Data Variability
  variability = check_variability(data = subdata)
  if (variability[['flag']]){
    return(variability)
    stop('Zero-variability variable detected, see variability object for details.')
  }
  
  # Calculate Data Thresholds
  
  if (is.null(penalties)){
    initial_penalties = array(data = rep(1, numvars*numvars), 
                              dim = c(numvars, numvars),
                              dimnames = list(c(colnames(subdata[[1]])),
                                              c(colnames(subdata[[1]]))))
  } else {
    initial_penalties = penalties
    colnames(initial_penalties) = colnames(subdata[[1]])
    rownames(initial_penalties) = colnames(subdata[[1]])
  }
  
  # Free AR Paths if Desired
  for (varname in yvarnames){
    if (ar == TRUE){
      initial_penalties[paste0(varname,'Lag'), varname] = 0 
    }
  }
  
  # Return Sample Penalty Matrix and Exit if Needed
  if (test_penalties == TRUE){
    return(initial_penalties)
    stop('Returning sample penatly matrix and exiting.')
  }
  
  # Group level search 
  grppaths <- group_search(subdata,
                           groupcutoff,
                           varname,
                           interact_exogenous,
                           predict_with_interactions,
                           interactnames,
                           interact_exogvars, 
                           output,
                           grppen = NULL)
  
  # Loop Through Subjects Again with the Group Level Information
  finalpaths <- ind_search(subdata,
                           varname,
                           interact_exogenous,
                           predict_with_interactions,
                           interactnames,
                           interact_exogvars, 
                           grppen = grppaths$group_penalties)
  
  # Optional search for subgroups using results from above. 
  if(subgroup){
    subgroup_results <- subgroup_search(subdata, 
                                        indpaths = finalpaths, 
                                        sub_method)
    if(subgroup_results$n_subgroups>1){
      subgrouppaths <- list()
      indpaths_sub  <- list()
      for (j in 1:subgroup_results$n_subgroups){
        sub_s_subjids = subset(subgroup_results$sub_mem$names,
                               subgroup_results$sub_mem$sub_membership == j)
        subgroupdata = subdata[sub_s_subjids]
        
        subgrouppaths[[j]] <- group_search(subdata = subgroupdata,
                                           groupcutoff,
                                           varname,
                                           interact_exogenous,
                                           predict_with_interactions,
                                           interactnames,
                                           interact_exogvars, 
                                           output, 
                                           grppen = grppaths$group_penalties)
        
        indpaths_sub[[j]] <- ind_search(subdata = subgroupdata,
                                        varname,
                                        interact_exogenous,
                                        predict_with_interactions,
                                        interactnames,
                                        interact_exogvars, 
                                        grppen = subgrouppaths[[j]]$group_penalties)  
        # combine subgroup-specific resuts into final estimates
        for (sub in names(subgroupdata)){
          finalpaths[,,sub] <- indpaths_sub[[j]][,,sub]
        }
      }
    }
  }

  
  # Organize Output
  grppaths$group_thresh_mat[is.na(grppaths$group_thresh_mat)] = 0
  output[['group']][['group_paths_present']] = grppaths$group_thresh_mat
  output[['group']][['group_penalties']] = grppaths$group_penalties
  for (sub in names(subdata)){
    output[[sub]][['data']] = subdata[[sub]]
    output[[sub]][['regression_matrix']] = finalpaths[, , sub]
  }
  
  if(subgroup)
    for (j in 1:subgroup_results$n_subgroups){
      output[['subgroup']][['subgroup_paths_present']][[j]] = subgrouppaths[[j]]$group_thresh_mat
      output[['subgroup']][['subgroup_group_penalties']][[j]] = subgrouppaths[[j]]$group_penalties
    } 
  
  # Add Visualization
  if (plot){
    output = network_vis(output)
  }
  
  # Save Output to Files
  manage_output(out = out, plot = plot, output = output)
  
  print('Algorithm successfully completed.', quote = FALSE)
  return(output)
}
