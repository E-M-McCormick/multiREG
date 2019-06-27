#' @name network_visualization
#' @aliases network_vis network_visualization
#' @title Network Visualization following Regularized Search
#' @description This function utilizes output from network_reg to build network plots.
#' @usage 
#' network_vis(output = NULL)
#'
#'@param output Output object from network_reg
network_visualization = network_vis = function(output = NULL){
  library(qgraph)
  
  # Create Individual-Level Figures
  subnames = names(output[!names(output) %in% c('group','function_parameters','variablenames')])
  for (sub in subnames){
    temp_regmat = output[[sub]][['regression_matrix']]
    contemp = temp_regmat[!grepl('Lag', rownames(temp_regmat)) & !grepl('_by_', rownames(temp_regmat)),]
    lagged = temp_regmat[grepl('Lag', rownames(temp_regmat)) & !grepl('_by_', rownames(temp_regmat)),]
    interacted = temp_regmat[grepl('_by_', rownames(temp_regmat)),]
    
    con.list = cbind(which(contemp!=0 & !is.na(contemp), arr.ind=TRUE), contemp[contemp!=0 & !is.na(contemp)])
    lag.list = cbind(which(lagged!=0 & !is.na(lagged), arr.ind=TRUE), lagged[lagged!=0 & !is.na(lagged)])
    int.list = cbind(which(interacted!=0 & !is.na(interacted), arr.ind=TRUE), interacted[interacted!=0 & !is.na(interacted)])
    edge.list = rbind(con.list, lag.list)
    
    is_exogenous = rownames(contemp) %in% output[['variablenames']][['exogenous_vars']]
    is_lagged = grepl('Lag',rownames(edge.list)) & !grepl('_by_', rownames(edge.list))
    rownames(edge.list)[is_lagged] = sub('Lag','',rownames(edge.list)[is_lagged])
    output[[sub]][['main_effects_fig']]= qgraph(edge.list,
                                                layout='circle',
                                                lty=ifelse(is_lagged,2,1),
                                                posCol='red',
                                                negCol='blue',
                                                esize = 5,
                                                parallelEdge = TRUE,
                                                fade = FALSE,
                                                labels = sub('_',' ', rownames(contemp)),
                                                label.cex = 1,
                                                shape = ifelse(is_exogenous,'square','circle'),
                                                DoNotPlot=TRUE)
    
    if (nrow(int.list) > 0){
      int.list = cbind(int.list,seq.int(nrow(int.list)))
      int.list2 = int.list1 = int.list
      rownames(int.list1) = sub('_by_.*','',rownames(int.list))
      rownames(int.list2) = sub('.*_by_','',rownames(int.list))
      
      is_lagged1 = grepl('Lag',rownames(int.list1))
      is_lagged2 = grepl('Lag',rownames(int.list2))
      
      rownames(int.list1)[is_lagged1] = sub('Lag','',rownames(int.list1)[is_lagged1])
      rownames(int.list2)[is_lagged2] = sub('Lag','',rownames(int.list2)[is_lagged2])
      
      int.list1[,1] = which(!is.na(contemp[,1]), arr.ind=T)[rownames(int.list1)]
      int.list2[,1] = which(!is.na(contemp[,1]), arr.ind=T)[rownames(int.list2)]
      
      int.list.final = rbind(int.list1, int.list2)
      is_lagged = c(is_lagged1, is_lagged2)
      
      output[[sub]][['interaction_fig']] = qgraph(int.list.final[,1:3],
                                                  layout='circle',
                                                  lty=ifelse(is_lagged,2,1),
                                                  posCol='red',
                                                  negCol='blue',
                                                  esize = 5,
                                                  parallelEdge = TRUE,
                                                  fade = FALSE,
                                                  labels = sub('_',' ', rownames(contemp)),
                                                  label.cex = 1,
                                                  shape = ifelse(is_exogenous,'square','circle'),
                                                  knots = int.list.final[,4],
                                                  knot.size = 3,
                                                  knot.color = 'green',
                                                  knot.borders = TRUE,
                                                  knot.border.width = 1,
                                                  DoNotPlot=TRUE)
    } else {
      output[[sub]][['interaction_fig']] = qgraph(int.list,
                                                  layout='circle',
                                                  posCol='red',
                                                  negCol='blue',
                                                  esize = 5,
                                                  parallelEdge = TRUE,
                                                  fade = FALSE,
                                                  labels = sub('_',' ', rownames(contemp)),
                                                  label.cex = 1,
                                                  shape = ifelse(is_exogenous,'square','circle'),
                                                  DoNotPlot=TRUE)
    }
  }
  # Create Subgroup-Level Figures
  
  
  # Create Group-Level Figure
  groupcut = output[["function_parameters"]][["groupcutoff"]]
  temp_counts = output[["group"]][["group_paths_proportions"]]
  contemp = temp_counts[!grepl('Lag', rownames(temp_counts)) & !grepl('_by_', rownames(temp_counts)),]
  lagged = temp_counts[grepl('Lag', rownames(temp_counts)) & !grepl('_by_', rownames(temp_counts)),]
  interacted = temp_counts[grepl('_by_', rownames(temp_counts)),]
  present = output[["group"]][["group_paths_present"]][!grepl('Lag', rownames(temp_counts)) & !grepl('_by_', rownames(temp_counts)),]
  
  con.list = cbind(which(contemp!=0 & !is.na(contemp), arr.ind=TRUE), contemp[contemp!=0 & !is.na(contemp)])
  lag.list = cbind(which(lagged!=0 & !is.na(lagged), arr.ind=TRUE), lagged[lagged!=0 & !is.na(lagged)])
  int.list = cbind(which(interacted!=0 & !is.na(interacted), arr.ind=TRUE), interacted[interacted!=0 & !is.na(interacted)])
  edge.list = rbind(con.list, lag.list)
  
  is_exogenous = rownames(contemp) %in% output[['variablenames']][['exogenous_vars']]
  is_lagged = grepl('Lag',rownames(edge.list)) & !grepl('_by_', rownames(edge.list))
  is_group = edge.list[,3] >= groupcut
  #is_subgroup = TO BE ADDED
  rownames(edge.list)[is_lagged] = sub('Lag','',rownames(edge.list)[is_lagged])

  output[['group']][['main_effects_fig']]= qgraph(edge.list,
                                              layout='circle',
                                              lty=ifelse(is_lagged,2,1),
                                              esize = 5,
                                              parallelEdge = TRUE,
                                              fade = FALSE,
                                              labels = sub('_',' ', rownames(contemp)),
                                              label.cex = 1,
                                              shape = ifelse(is_exogenous,'square','circle'),
                                              edge.width = ifelse(is_group, 1, edge.list[,3]),
                                              #edge.width = ifelse(is_group,1,ifelse(is_subgroup,1,edge.list[,3])),
                                              edge.color = ifelse(is_group,'black','grey'),
                                              #edge.color = ifelse(is_group,'black',ifelse(is_subgroup,'green','grey')),
                                              DoNotPlot=TRUE)
  if (nrow(int.list) > 0){
    int.list = cbind(int.list,seq.int(nrow(int.list)))
    int.list2 = int.list1 = int.list
    rownames(int.list1) = sub('_by_.*','',rownames(int.list))
    rownames(int.list2) = sub('.*_by_','',rownames(int.list))
    
    is_lagged1 = grepl('Lag',rownames(int.list1))
    is_lagged2 = grepl('Lag',rownames(int.list2))
    
    rownames(int.list1)[is_lagged1] = sub('Lag','',rownames(int.list1)[is_lagged1])
    rownames(int.list2)[is_lagged2] = sub('Lag','',rownames(int.list2)[is_lagged2])
    
    int.list1[,1] = which(!is.na(present[,1]), arr.ind=T)[rownames(int.list1)]
    int.list2[,1] = which(!is.na(present[,1]), arr.ind=T)[rownames(int.list2)]
    
    int.list.final = rbind(int.list1, int.list2)
    is_lagged = c(is_lagged1, is_lagged2)
    is_group = int.list.final[,3] >= groupcut
    
    output[['group']][['interaction_fig']] = qgraph(int.list.final[,1:3],
                                                layout='circle',
                                                lty=ifelse(is_lagged,2,1),
                                                esize = 5,
                                                parallelEdge = TRUE,
                                                fade = FALSE,
                                                labels = sub('_',' ', rownames(contemp)),
                                                label.cex = 1,
                                                shape = ifelse(is_exogenous,'square','circle'),
                                                edge.width = ifelse(is_group, 1, int.list.final[,3]),
                                                #edge.width = ifelse(is_group,1,ifelse(is_subgroup,1,edge.list[,3])),
                                                edge.color = ifelse(is_group, 'black', 'grey'),
                                                #edge.color = ifelse(is_group,'black',ifelse(is_subgroup,'green','grey')),
                                                knots = int.list.final[,4],
                                                knot.size = ifelse(is_group,3,1),
                                                knot.color = ifelse(is_group, 'green', 'grey'),
                                                knot.borders = TRUE,
                                                knot.border.width = 1,
                                                DoNotPlot=TRUE)
    
  } else {
    output[['group']][['interaction_fig']] = qgraph(int.list,
                                                layout='circle',
                                                posCol='red',
                                                negCol='blue',
                                                esize = 5,
                                                parallelEdge = TRUE,
                                                fade = FALSE,
                                                labels = sub('_',' ', rownames(contemp)),
                                                label.cex = 1,
                                                shape = ifelse(is_exogenous,'square','circle'),
                                                DoNotPlot=TRUE)
  }
  return(output)
}
