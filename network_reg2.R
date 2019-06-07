network_reg <- function(data = '',
                        sep = '',
                        header = TRUE,
                        group_cutoff = .75,
                        out = '',
                        alpha = .5,
                        ar = TRUE,
                        penalties = NULL,
                        exogenous = NULL,
                        lag_exogenous = FALSE){
  
  library(tools); library(glmnet) 
  refpath = getwd()
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
  
  # Categorize Variables. & Omit NaN Rows
  for (i in 1:length(subdata)){
    yvar = subdata[[i]][,!colnames(subdata[[i]]) %in% exogenous, drop=FALSE]
    yvarnames = colnames(yvar)
    lagvar = rbind(rep(NA, ncol(yvar)), yvar[1:(nrow(yvar)-1),])
    colnames(lagvar) = paste(colnames(lagvar), 'Lag', sep='')
    exogvar = subdata[[i]][,colnames(subdata[[i]]) %in% exogenous, drop=FALSE]
    if (lag_exogenous == TRUE){
      lagexogvar = rbind(rep(NA, ncol(exogvar)), exogvar[1:(nrow(exogvar)-1), , drop=FALSE])
      colnames(lagexogvar) = paste(colnames(exogvar), 'Lag', sep='')
      subdata[[i]] = na.omit(cbind(yvar,lagvar,exogvar,lagexogvar))
    } else if (lag_exogenous == FALSE){
      subdata[[i]] = na.omit(cbind(yvar,lagvar,exogvar))
    }
  }
  print('Data successfully read in.')
  
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
  initial_penalties = array(data = rep(1, numvars*numvars), 
                            dim = c(numvars,numvars),
                            dimnames = list(c(colnames(subdata[[1]])),
                                            c(colnames(subdata[[1]]))))
  # Free AR Paths if Desired
  for (varname in yvarnames){
    if (ar == TRUE){
      initial_penalties[paste0(varname,'Lag'),varname] = 0 
    }
  }
  
  # Loop through Subjects Data to Build Individual Models
  for (sub in names(subdata)){
    print(paste0('Building group-level model for ', sub, '.'))
    tempdata = subdata[[sub]]
    for (varname in yvarnames){
      cvfit = cv.glmnet(x = as.matrix(tempdata[,!colnames(tempdata) %in% varname]),
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
  output = list()
  group_thresh_mat = rowSums(pathpresent, dims = 2)
  output[['group']][['group_paths_counts']] = group_thresh_mat
  
  group_thresh_mat = group_thresh_mat/nsubs
  group_thresh_mat[group_thresh_mat < group_cutoff] = 0
  group_thresh_mat[group_thresh_mat >= group_cutoff] = 1
  group_penalties = abs(group_thresh_mat - 1)
  
  # Loop Through Subjects Again with the Group Level Information
  for (sub in names(subdata)){
    print(paste0('Building individual-level model for ', sub, '.'))
    tempdata = subdata[[sub]]
    for (varname in yvarnames){
      cvfit = cv.glmnet(x = as.matrix(tempdata[,!colnames(tempdata) %in% varname]),
                        y = tempdata[,colnames(tempdata) %in% varname],
                        type.measure = 'mse',
                        nfolds = round(nrow(tempdata)/10),
                        alpha = alpha,
                        penalty.factor = group_penalties[!colnames(tempdata) %in% varname, varname])
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
  
  # Save Output to Files
  if (out != ''){
    print('Writing output to file.')
    if (!dir.exists(out)){
      print('Creating output directories')
      dir.create(out);
      dir.create(paste(out, 'individual', sep=.Platform$file.sep))
    }
    write.csv(output$group$group_paths_counts[, colnames(output$group$group_paths_counts) %in% yvarnames],
              file = paste(out, 'GroupLevel_PathCountsMatrix.csv', sep=.Platform$file.sep))
    write.csv(output$group$group_paths_present[, colnames(output$group$group_paths_present) %in% yvarnames],
              file = paste(out, 'GroupLevel_PathsPresent.csv', sep=.Platform$file.sep))
    for (sub in names(subdata)){
      write.csv(output[[sub]]$regression_matrix[, colnames(output$group$group_paths_counts) %in% yvarnames],
                file = paste(out, 'individual', paste0(sub,'_Betas.csv'), sep=.Platform$file.sep))
    }
  }
  setwd(refpath)
  return(output)
}