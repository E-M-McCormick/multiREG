#' @name netreg
#' @aliases network_reg netreg
#' @title Network Regularization Model Search
#' @description This function utilizes regression with regularization to build models for individuals 
#' consisting of individual and group-level paths.
#' @usage 
#' network_reg(data                    = NULL,
#'             out                     = NULL,
#'             sep                     = NULL,
#'             header                  = TRUE,
#'             ar                      = TRUE,
#'             plot                    = TRUE,
#'             subgroup                = FALSE,
#'             sub_feature             = 'lag & contemp',       
#'             sub_method              = 'walktrap',
#'             confirm_subgroup        = NULL,
#'             paths                   = NULL,
#'             conv_vars               = NULL,
#'             conv_length             = 16,
#'             conv_interval           = 1,
#'             mult_vars               = NULL,
#'             mean_center_mult        = FALSE,
#'             standardize             = FALSE,
#'             groupcutoff             = .75,
#'             subcutoff               = .5,
#'             diagnos                 = FALSE,
#'             alpha                   = .5,
#'             penalties               = NULL,
#'             test_penalties          = FALSE,
#'             exogenous               = NULL,
#'             lag_exogenous           = FALSE,
#'             interact_exogenous      = NULL,
#'             interact_with_exogenous = NULL,
#'             predict_with_interactions  = NULL)
#'             
#' @param data The path to the directory where individual data files are located,
#' or the name of the list containing individual data. Each file or matrix within the list
#' must contain a single matrix containing the a T (time) by p (number of variables) matrix,
#' where the rows represent time and columns represent individual variables. Individuals may
#' have different numbers of observations (T), but must have the same number of variables (p).
#' 
#' @param sep Spacing scheme for input files. 
#' '' indicates space-separated; ',' indicates comma separated; '/t' indicates tab-separated
#' Only necessary when reading in files from physical directory.
#' 
#' @param header (Logical) Indicate TRUE if variable names incluced in input file, FALSE otherwise.
#' Only necessary when reading in files from physical directory.
#' 
#' @param ar (Logical) If TRUE, begin model search with all autoregressive pathways estimated
#' with no shrinkage (i.e., penalty = 0).
#' 
#' @param group_cutoff Cutoff value for inclusion of a given path at the group-level.
#' For instance, group_cutoff = .75 indicates that a path needs to be estimated for 75% of
#' individuals to be included as a group-level path.
#' 
#' @param out (Optional) The path to directory where results will be stored. If specified,
#' a copy of output data will be saved into the directory. If the specified directory does
#' not exist, it will be created.
#' 
#' @param alpha Elastic-net parameter for the regularization approach. Values close to 0 mimic 
#' the ridge penalty, which tends to shrink correlated parameters towards one another. Values 
#' close to 1 mimic the lasso penalty, which tends to select one parameter and shrink
#' the others. The default value (alpha=.5) balances these two considerations, and tends to select
#' groups of correlated parameters and shrink other groups towards zero.
#' 
#' @param penalties (Optional) A matrix of user-provided penalties to initialize group-model search. 
#' Should contain a column for all variables (including lagged versions and interactions) that will 
#' be included in the model search. Values of 1 (the default) will initilize a variable to be 
#' normally considered in the regularization, values of 0 will initilize a variable to be estimated
#' (i.e., no shrinkage), and values of Inf will exlcude variables from the model.
#' 
#'  @param test_penalties (Optional, Logical) Optional argument to output a sample penalty matrix
#'  based on function parameters. Helpful for specifying a matrix to use in the penalties arguement.
#'  Function will exit gracefully before running anything if test_penalties = TRUE.
#' 
#' @param exogenous (Optional) A list of user-specified variables to consider as exogenous
#' (e.g., cannot be predicted) in the model search procedure. If variable names are supplied,
#' variables should be referred to by name. If not, then variables should be referenced by
#' the pattern 'V#', where # represents the column number in the original data file (e.g., 'V5').
#' 
#' @param lag_exogenous (Optional, Logical) If TRUE, a lagged version of the exogenous variable 
#' will be created. If set to TRUE, but exogenous variables are not indicated in the arguement
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

network_reg <- netreg <- function(data                    = NULL,
                                  out                     = NULL,
                                  sep                     = NULL,
                                  header                  = TRUE,
                                  ar                      = TRUE,
                                  plot                    = TRUE,
                                  # subgroup                = FALSE,                 I added these arguements to be consistent with the gimme code, but so far they are not supported.
                                  # sub_feature             = 'lag & contemp',       
                                  # sub_method              = 'walktrap',
                                  # confirm_subgroup        = NULL,
                                  # paths                   = NULL,
                                  conv_vars               = NULL,
                                  conv_length             = 16,
                                  conv_interval           = 1,
                                  # mult_vars               = NULL,
                                  # mean_center_mult        = FALSE,
                                  # standardize             = FALSE,
                                  groupcutoff             = .75,
                                  # subcutoff               = .5,
                                  # diagnos                 = FALSE,                Might not include depending on what this does. Couldn't tell from the code I looked at.
                                  alpha                   = .5,
                                  penalties               = NULL,
                                  test_penalties          = FALSE,
                                  exogenous               = NULL,
                                  lag_exogenous           = FALSE,
                                  interact_exogenous      = NULL,
                                  interact_with_exogenous = NULL,
                                  predict_with_interactions  = NULL){
  
  
  library(tools); library(glmnet); library(gimme)

  # Add Function Parameters to Output
  output = list()
  output[['function_parameters']] = as.list(sys.call())
  
  # Wrangle Data into List
  if (!is.list(data)){
    subdata = list(); setwd(data)
    for (i in list.files(data)){
      tempname = file_path_sans_ext(i)
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
  interact_with_all_flag = FALSE
  if (!is.null(interact_with_exogenous) && interact_with_exogenous == 'all'){
      interact_with_all_flag = TRUE
  }
  for (i in 1:length(subdata)){
    yvar = subdata[[i]][,!colnames(subdata[[i]]) %in% exogenous, drop=FALSE]
    yvarnames = colnames(yvar)
    lagvar = rbind(rep(NA, ncol(yvar)), yvar[1:(nrow(yvar)-1),])
    colnames(lagvar) = paste(colnames(lagvar), 'Lag', sep='')
    exogvar = subdata[[i]][,colnames(subdata[[i]]) %in% exogenous, drop=FALSE]
    
    
    # Lag Exogenous Variables if Needed
    if (lag_exogenous == TRUE){
      lagexogvar = rbind(rep(NA, ncol(exogvar)), exogvar[1:(nrow(exogvar)-1), , drop=FALSE])
      colnames(lagexogvar) = paste(colnames(exogvar), 'Lag', sep='')
      if (!is.null(interact_exogenous) && interact_exogenous == 'all'){
        interact_exogvars = c(colnames(exogvar),colnames(lagexogvar))
      } else {
        interact_exogvars = interact_exogenous
      }
      if (!is.null(interact_exogenous)){
        
        # Create Interaction Variables As Needed
        if (interact_with_all_flag == TRUE){
          interact_var = sapply(cbind(yvar,lagvar), 
                                function(x){ x*cbind(exogvar[colnames(exogvar) %in% interact_exogvars],lagexogvar[colnames(lagexogvar) %in% interact_exogvars]) })
          newnames = as.vector(sapply(cbind(colnames(yvar),colnames(lagvar)), 
                                      function(x){ paste0(x,'_x_',cbind(colnames(exogvar[colnames(exogvar) %in% interact_exogvars]),colnames(lagexogvar[colnames(lagexogvar) %in% interact_exogvars]))) }))
        } else {
          
          #####
          interact_var = sapply(cbind(yvar[colnames(yvar) %in% interact_with_exogenous],lagvar[colnames(lagvar) %in% interact_with_exogenous]), 
                                function(x){ x*cbind(exogvar[colnames(exogvar) %in% interact_exogvars],lagexogvar[colnames(lagexogvar) %in% interact_exogvars]) })
          
          
          newnames = unique(as.vector(sapply(cbind(colnames(yvar[colnames(yvar) %in% interact_with_exogenous]),colnames(lagvar[colnames(lagvar) %in% interact_with_exogenous])), 
                                             function(x){ paste0(x,'_x_',cbind(colnames(exogvar[colnames(exogvar) %in% interact_exogvars]),colnames(lagexogvar[colnames(lagexogvar) %in% interact_exogvars]))) })))
          
          #####
        }
        interact_var = matrix(unlist(interact_var), ncol=length(interact_var))
        colnames(interact_var) = newnames
        interactnames = newnames
        subdata[[i]] = na.omit(cbind(yvar,lagvar,exogvar,lagexogvar,interact_var))
      } else {
        interact_exogvars = ''
        interactnames = ''
        subdata[[i]] = na.omit(cbind(yvar,lagvar,exogvar,lagexogvar))
      }
      exognames = c(colnames(exogvar),colnames(lagexogvar))
    } else if (lag_exogenous == FALSE){
      if (!is.null(interact_exogenous) && interact_exogenous == 'all'){
        interact_exogvars = colnames(exogvar)
      } else {
        interact_exogvars = interact_exogenous
      }
      # Create Interaction Variables As Needed
      if (!is.null(interact_exogenous)){
        if (interact_with_all_flag == TRUE){
          interact_var = sapply(cbind(yvar,lagvar), 
                                function(x,y){ x*exogvar[colnames(exogvar) %in% interact_exogvars] })
          newnames = as.vector(sapply(cbind(colnames(yvar),colnames(lagvar)), 
                                      function(x){ paste0(x,'_x_',colnames(exogvar[colnames(exogvar) %in% interact_exogvars])) }))
        } else {
          interact_var = sapply(cbind(yvar[colnames(yvar) %in% interact_with_exogenous],lagvar[colnames(lagvar) %in% interact_with_exogenous]), 
                                function(x){ x*exogvar[colnames(exogvar) %in% interact_exogvars] })
          newnames = unique(as.vector(sapply(cbind(colnames(yvar[colnames(yvar) %in% interact_with_exogenous]),colnames(lagvar[colnames(lagvar) %in% interact_with_exogenous])), 
                                             function(x){ paste0(x,'_x_',colnames(exogvar[colnames(exogvar) %in% interact_exogvars])) })))
        }
        interact_var = matrix(unlist(interact_var), ncol=length(interact_var))
        colnames(interact_var) = newnames
        interactnames = newnames
        subdata[[i]] = na.omit(cbind(yvar,lagvar,exogvar,interact_var))
      } else {
        interact_exogvars = ''
        interactnames = ''
        subdata[[i]] = na.omit(cbind(yvar,lagvar,exogvar))
      }
      exognames = colnames(exogvar)
    }
  }
  print('Data successfully read in.')
  output[['variablenames']][['y_vars']] = yvarnames
  output[['variablenames']][['exogenous_vars']] = exognames
  output[['variablenames']][['interaction_vars']] = interactnames
  
  # Calculate Data Thresholds
  nsubs = length(subdata)
  numvars = ncol(subdata[[1]])
  pathpresent = array(data = rep(NaN, numvars*numvars*length(subdata)), 
                      dim = c(numvars,numvars,length(subdata)),
                      dimnames = list(c(colnames(subdata[[1]])),
                                      c(colnames(subdata[[1]])),
                                      c(names(subdata))))
  finalpaths = array(data = rep(0, numvars*numvars*length(subdata)), 
                     dim = c(numvars,numvars,length(subdata)),
                     dimnames = list(c(colnames(subdata[[1]])),
                                     c(colnames(subdata[[1]])),
                                     c(names(subdata))))
  if (is.null(penalties)){
    initial_penalties = array(data = rep(1, numvars*numvars), 
                              dim = c(numvars,numvars),
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
      initial_penalties[paste0(varname,'Lag'),varname] = 0 
    }
  }
  
  # Return Sample Penalty Matrix and Exit if Needed
  if (test_penalties == TRUE){
    return(initial_penalties)
    stop('Returning sample penatly matrix and exiting.')
  }
  
  # Loop through Subjects Data to Build Individual Models
  for (sub in names(subdata)){
    print(paste0('Building group-level model for ', sub, '.'))
    tempdata = subdata[[sub]]
    for (varname in yvarnames){
      subset_predictors = as.matrix(tempdata[,!(colnames(tempdata) %in% varname |
                                                  colnames(tempdata) %in% paste0(varname,'_x_',interact_exogvars))])
      if (!is.null(predict_with_interactions) & !varname %in% predict_with_interactions){
          subset_predictors = subset_predictors[,!colnames(subset_predictors) %in% interactnames]
      }
      cvfit = cv.glmnet(x = subset_predictors,
                        y =  tempdata[,colnames(tempdata) %in% varname],
                        type.measure = 'mse',
                        nfolds = round(nrow(tempdata)/10),
                        alpha = alpha,
                        penalty.factor = initial_penalties[!colnames(tempdata) %in% varname, varname])
      temp = coef(cvfit, s = 'lambda.1se')
      for (predictor in rownames(temp)[!rownames(temp) %in% '(Intercept)']){
        if (temp[predictor,] == 0){
          pathpresent[predictor,varname,sub] = 0
        } else {
          pathpresent[predictor,varname,sub] = 1
        }
      }
    }
  }
  
  # Calculate Paths that Should Appear in the Group (Non-Penalized) Model
  group_thresh_mat = rowSums(pathpresent, dims = 2)
  output[['group']][['group_paths_counts']] = group_thresh_mat
  
  group_thresh_mat = group_thresh_mat/nsubs
  output[['group']][['group_paths_proportions']] = group_thresh_mat
  group_thresh_mat[group_thresh_mat < groupcutoff] = 0
  group_thresh_mat[group_thresh_mat >= groupcutoff] = 1
  group_penalties = abs(group_thresh_mat - 1)
  
  # Loop Through Subjects Again with the Group Level Information
  for (sub in names(subdata)){
    print(paste0('Building individual-level model for ', sub, '.'))
    tempdata = subdata[[sub]]
    for (varname in yvarnames){
      subset_predictors = as.matrix(tempdata[,!(colnames(tempdata) %in% varname |
                                                  colnames(tempdata) %in% paste0(varname,'_x_',interact_exogvars))])
      if (!is.null(predict_with_interactions) & !varname %in% predict_with_interactions){
        subset_predictors = subset_predictors[,!colnames(subset_predictors) %in% interactnames]
      }
      subset_penalties = group_penalties[!(colnames(tempdata) %in% varname |
                                              colnames(tempdata) %in% paste0(varname,'_x_',interact_exogvars)), varname]
      if (!is.null(predict_with_interactions) & !varname %in% predict_with_interactions){
        subset_penalties = subset_penalties[!names(subset_penalties) %in% interactnames]
      }
      cvfit = cv.glmnet(x = subset_predictors,
                        y = tempdata[,colnames(tempdata) %in% varname],
                        type.measure = 'mse',
                        nfolds = round(nrow(tempdata)/10),
                        alpha = alpha,
                        penalty.factor = subset_penalties)
      temp = coef(cvfit, s = 'lambda.1se')
      for (predictor in rownames(temp)[!rownames(temp) %in% '(Intercept)']){
        if (temp[predictor,] != 0){
          finalpaths[predictor,varname,sub] = temp[predictor,]
        }
      }
    }
  }
  
  # Organize Output
  group_thresh_mat[is.na(group_thresh_mat)] = 0
  output[['group']][['group_paths_present']] = group_thresh_mat
  output[['group']][['group_penalties']] = group_penalties
  for (sub in names(subdata)){
    output[[sub]][['data']] = subdata[[sub]]
    output[[sub]][['regression_matrix']] = finalpaths[, , sub]
  }
  
  # Add Visualization
  if (plot){
    output = network_vis(output)
  }
  
  
  # Save Output to Files
  if (!is.null(out)){
    if (!dir.exists(out)){
      print('Creating output directories')
      dir.create(out);
    }
    print('Writing output to file.')
    write.csv(output$group$group_paths_counts[, colnames(output$group$group_paths_counts) %in% yvarnames],
              file = paste(out, 'GroupLevel_PathCountsMatrix.csv', sep=.Platform$file.sep))
    write.csv(output$group$group_paths_present[, colnames(output$group$group_paths_present) %in% yvarnames],
              file = paste(out, 'GroupLevel_PathsPresent.csv', sep=.Platform$file.sep))
    if (plot) {
      pdf(file.path(out, 'GroupLevel_main_effects_plot.pdf'))
      plot(output[['group']][['main_effects_fig']])
      dev.off()
      pdf(file.path(out, 'GroupLevel_interactions_plot.pdf'))
      plot(output[['group']][['interaction_fig']])
      dev.off()
    }
    
    
    dir.create(paste(out, 'individual', sep=.Platform$file.sep))
    for (sub in names(subdata)){
      write.csv(output[[sub]]$regression_matrix[, colnames(output$group$group_paths_counts) %in% yvarnames],
                file = paste(out, 'individual', paste0(sub,'_Betas.csv'), sep=.Platform$file.sep))
      if (plot) {
        pdf(file.path(out, 'individual', paste0(sub,'_main_effects_plot.pdf')))
        plot(output[[sub]][['main_effects_fig']])
        dev.off()
        pdf(file.path(out, 'individual', paste0(sub,'_interactions_plot.pdf')))
        plot(output[[sub]][['interaction_fig']])
        dev.off()
      }
    }
    setwd(out)
  }
  return(output)
}
