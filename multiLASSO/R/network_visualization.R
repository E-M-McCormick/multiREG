#' @name network_visualization
#' @aliases network_vis network_visualization
#' @title Network Visualization Following Regularization with LASSO
#' @keywords internal 
network_visualization = network_vis = function(output = NULL){
  
  # Determine if Interaction Plots are Needed
  interaction_flag = TRUE
  if (is.null(output[['function_parameters']][['interact_exogenous']]) & 
      is.null(output[['function_parameters']][['interact_with_exogenous']])) {
    interaction_flag = FALSE
  }
  
  # Create Individual-Level Figures
  subnames = names(output[!names(output) %in% c('group','function_parameters','variablenames')])
  if (length(output[['variablenames']][['y_vars']]) <= 100) {
    for (sub in subnames){
      print(paste0('Creating plots for ', sub, '.'), quote = FALSE)
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
      
      output[[sub]][['main_effects_fig']] = qgraph::qgraph(edge.list,
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
                                                           DoNotPlot=FALSE)
      
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
        
        output[[sub]][['interaction_fig']] = qgraph::qgraph(int.list.final[,1:3],
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
                                                            DoNotPlot=FALSE)
      } else {
        if (interaction_flag) {
          output[[sub]][['interaction_fig']] = qgraph::qgraph(int.list,
                                                              layout='circle',
                                                              posCol='red',
                                                              negCol='blue',
                                                              esize = 5,
                                                              parallelEdge = TRUE,
                                                              fade = FALSE,
                                                              labels = sub('_',' ', rownames(contemp)),
                                                              label.cex = 1,
                                                              shape = ifelse(is_exogenous,'square','circle'),
                                                              DoNotPlot=FALSE)
        }
      }
    }
  } else {
    print('Suppressing individual-level plots to improve clarity and runtime.', quote = FALSE)
  }
  # Create Subgroup-Level Figures
  
  
  # Create Group-Level Figure
  print('Creating Group Plots.', quote = FALSE)
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

  if (length(output[['variablenames']][['y_vars']]) <= 100) {
    output[['group']][['main_effects_fig']] = qgraph::qgraph(edge.list,
                                                             layout='circular',
                                                             lty=ifelse(is_lagged,2,1),
                                                             esize = 3,
                                                             parallelEdge = TRUE,
                                                             fade = FALSE,
                                                             labels = sub('_',' ', rownames(contemp)),
                                                             label.cex = 1,
                                                             shape = ifelse(is_exogenous,'square','circle'),
                                                             edge.width = ifelse(is_group, 1, edge.list[,3]),
                                                             #edge.width = ifelse(is_group,1,ifelse(is_subgroup,1,edge.list[,3])),
                                                             edge.color = ifelse(is_group,'black','grey'),
                                                             #edge.color = ifelse(is_group,'black',ifelse(is_subgroup,'green','grey')),
                                                             DoNotPlot=FALSE)
  } else {
    edge.list = edge.list[is_group, ]
    is_lagged = is_lagged[is_group]
    output[['group']][['main_effects_fig']] = qgraph::qgraph(edge.list[, 1:2],
                                                             layout = 'circle',
                                                             lty=ifelse(is_lagged,2,1),
                                                             esize = 1,
                                                             parallelEdge = TRUE,
                                                             labels = sub('_',' ', rownames(contemp)),
                                                             label.cex = 1,
                                                             shape = ifelse(is_exogenous,'square','circle'),
                                                             edge.width = 1,
                                                             edge.color = 'black',
                                                             DoNotPlot = FALSE)
  }        
                                                           
  
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
    
    if (length(output[['variablenames']][['y_vars']]) <= 100) {
      output[['group']][['interaction_fig']] = qgraph::qgraph(int.list.final[, 1:3],
                                                              layout='circle',
                                                              lty=ifelse(is_lagged,2,1),
                                                              esize = 3,
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
                                                              DoNotPlot=FALSE)
    } else {
      int.list.final = int.list.final[is_group, ]
      is_lagged = is_lagged[is_group]
      output[['group']][['interaction_fig']] = qgraph::qgraph(int.list.final[, 1:2],
                                                              layout='circle',
                                                              lty=ifelse(is_lagged,2,1),
                                                              esize = 1,
                                                              parallelEdge = TRUE,
                                                              fade = FALSE,
                                                              labels = sub('_',' ', rownames(contemp)),
                                                              label.cex = 1,
                                                              shape = ifelse(is_exogenous,'square','circle'),
                                                              edge.width = 1,
                                                              #edge.width = ifelse(is_group,1,ifelse(is_subgroup,1,edge.list[,3])),
                                                              edge.color = 'black',
                                                              #edge.color = ifelse(is_group,'black',ifelse(is_subgroup,'green','grey')),
                                                              knots = int.list.final[,4],
                                                              knot.size = 1,
                                                              knot.color = 'green',
                                                              knot.borders = TRUE,
                                                              knot.border.width = .5,
                                                              DoNotPlot=FALSE)
    }
    
  } else {
    if (interaction_flag) {
      output[['group']][['interaction_fig']] = qgraph::qgraph(int.list,
                                                              layout='circle',
                                                              posCol='red',
                                                              negCol='blue',
                                                              esize = 5,
                                                              parallelEdge = TRUE,
                                                              fade = FALSE,
                                                              labels = sub('_',' ', rownames(contemp)),
                                                              label.cex = 1,
                                                              shape = ifelse(is_exogenous,'square','circle'),
                                                              DoNotPlot=FALSE)
    }
  }
  dev.off()
  return(output)
}
