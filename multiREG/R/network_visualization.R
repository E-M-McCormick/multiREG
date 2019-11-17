#' Network Visualization
#' @param output Output object.
#' @param finalpaths Matrix of final coefficients for each path.
#' @param verbose Logical. If TRUE, algorithm will print progress to console.
#' @return Returns network plots.
#' @keywords internal 
network_visualization = network_vis = function(output = NULL, finalpaths = NULL, verbose = TRUE){
  
  #### Determine if Interaction Plots are Needed ####
  interaction_flag = TRUE
  if (is.null(output[['function_parameters']][['interact_exogenous']]) & 
      is.null(output[['function_parameters']][['interact_with_exogenous']])) {
    interaction_flag = FALSE
  }
  
  #### Create Individual-Level Figures ####
  subnames = names(output[!names(output) %in% c('group','function_parameters','variablenames')])
  if (output[['function_parameters']][['subgroup']]){subnames = subnames[subnames != 'subgroup']}
  if (length(output[['variablenames']][['y_vars']]) <= 100) {
    for (sub in subnames){
      if(verbose){print(paste0('Creating plots for ', sub, '.'), quote = FALSE)}
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
    if(verbose){print('Suppressing individual-level plots to improve clarity and runtime.', quote = FALSE)}
  }
  id_subpaths_all = matrix(0,length(output$group$group_paths_present[,1]), length(output$group$group_paths_present[1,]))
  
  #### Create Subgroup Plots if Needed ####
  if(output$function_parameters$subgroup && output$subgroup$subgroup_number>1 && length(subnames)>output$subgroup$subgroup_number && output$function_parameters$heuristic == 'GIMME'){
    # Create Subgroup-Level Figures, and get some info to add colors to the group-level plot
    moderated_group_by_subgroup = array(rep(0, nrow(id_subpaths_all)*ncol(id_subpaths_all)*output$subgroup$subgroup_number),
                                        dim = c(nrow(id_subpaths_all), ncol(id_subpaths_all), output$subgroup$subgroup_number),
                                        dimnames = list(c(colnames(id_subpaths_all)),
                                                        c(colnames(id_subpaths_all)),
                                                        c(1:output$subgroup$subgroup_number)))
    for (p in 1:output$subgroup$subgroup_number){
      if(verbose){print(paste0('Creating plots for subgroup', p, '.'), quote = FALSE)}
      if (length(which(output$subgroup$membership[,2]==p))>1){
        id_subpaths = output$subgroup$subgroup_paths_present[[p]] - output$group$group_paths_present
        sub_temp_counts = id_subpaths + output$subgroup$subgroup_paths_proportions[[p]] # subgroup-level paths will be >1
        
        # Test if any group paths differ significantly in sign
        group_subgroupsign = group_subgroupsign(output$group$group_paths_present, finalpaths, output$subgroup$membership)
        
        sub_temp_counts = sub_temp_counts*group_subgroupsign
        moderated_group_by_subgroup[,,p] = sub_temp_counts
        
        sub_contemp = sub_temp_counts[!grepl('Lag', rownames(sub_temp_counts)) & !grepl('_by_', rownames(sub_temp_counts)),]
        sub_lagged = sub_temp_counts[grepl('Lag', rownames(sub_temp_counts)) & !grepl('_by_', rownames(sub_temp_counts)),]
        sub_interacted = sub_temp_counts[grepl('_by_', rownames(sub_temp_counts)),]
        
        sub_con.list = cbind(which(sub_contemp!=0 & !is.na(sub_contemp), arr.ind=TRUE), sub_contemp[sub_contemp!=0 & !is.na(sub_contemp)])
        sub_lag.list = cbind(which(sub_lagged!=0 & !is.na(sub_lagged), arr.ind=TRUE), sub_lagged[sub_lagged!=0 & !is.na(sub_lagged)])
        sub_int.list = cbind(which(sub_interacted!=0 & !is.na(sub_interacted), arr.ind=TRUE), sub_interacted[sub_interacted!=0 & !is.na(sub_interacted)])
        sub_edge.list = rbind(sub_con.list, sub_lag.list)
        
        sub_is_exogenous = rownames(sub_contemp) %in% output[['variablenames']][['exogenous_vars']]
        sub_is_lagged = grepl('Lag',rownames(sub_edge.list)) & !grepl('_by_', rownames(sub_edge.list))
        is_subgroup = sub_edge.list[,3] > 1
        sub_is_group = dplyr::between(sub_edge.list[,3], output$function_parameters$groupcutoff, 1)
        group_is_sub = sign(sub_edge.list[,3]) == -1
        sub_edge.list[,3] = abs(sub_edge.list[,3])
        
        rownames(sub_edge.list)[sub_is_lagged] = sub('Lag','',rownames(sub_edge.list)[sub_is_lagged])
        
        output[['subgroup']][['main_effects_fig']][[p]] = qgraph::qgraph(sub_edge.list,
                                                                         layout='circle',
                                                                         lty=ifelse(sub_is_lagged,2,1),
                                                                         esize = 3,
                                                                         parallelEdge = TRUE,
                                                                         fade = FALSE,
                                                                         labels = sub('_',' ', rownames(sub_contemp)),
                                                                         label.cex = 1,
                                                                         shape = ifelse(sub_is_exogenous,'square','circle'),
                                                                         edge.width = ifelse(sub_is_group,1,ifelse(is_subgroup,1,sub_edge.list[,3])),
                                                                         edge.color = ifelse(sub_is_group,ifelse(group_is_sub,'purple','black'),ifelse(is_subgroup,'green','grey')),
                                                                         DoNotPlot=FALSE)
        
        if (nrow(sub_int.list) > 0){
          sub_int.list = cbind(sub_int.list,seq.int(nrow(sub_int.list)))
          sub_int.list2 = sub_int.list1 = sub_int.list
          rownames(sub_int.list1) = sub('_by_.*','',rownames(sub_int.list))
          rownames(sub_int.list2) = sub('.*_by_','',rownames(sub_int.list))
          
          sub_is_lagged1 = grepl('Lag',rownames(sub_int.list1))
          sub_is_lagged2 = grepl('Lag',rownames(sub_int.list2))
          
          rownames(sub_int.list1)[sub_is_lagged1] = sub('Lag','',rownames(sub_int.list1)[sub_is_lagged1])
          rownames(sub_int.list2)[sub_is_lagged2] = sub('Lag','',rownames(sub_int.list2)[sub_is_lagged2])
          
          sub_int.list1[,1] = which(!is.na(sub_contemp[,1]), arr.ind=T)[rownames(sub_int.list1)]
          sub_int.list2[,1] = which(!is.na(sub_contemp[,1]), arr.ind=T)[rownames(sub_int.list2)]
          
          sub_int.list.final = rbind(sub_int.list1, sub_int.list2)
          sub_is_lagged = c(sub_is_lagged1, sub_is_lagged2)
          is_subgroup_int = sub_int.list.final[,3] > 1
          sub_is_group_int = dplyr::between(sub_int.list.final[,3], output$function_parameters$groupcutoff, 1)
          group_is_sub_int = sign(sub_int.list.final[,3]) == -1
          
          output[['subgroup']][['interaction_fig']][[p]] = qgraph::qgraph(sub_int.list.final[,1:3],
                                                                          layout='circle',
                                                                          lty=ifelse(sub_is_lagged,2,1),
                                                                          esize = 1,
                                                                          parallelEdge = TRUE,
                                                                          fade = FALSE,
                                                                          labels = sub('_',' ', rownames(sub_contemp)),
                                                                          label.cex = 1,
                                                                          shape = ifelse(sub_is_exogenous,'square','circle'),
                                                                          edge.width = ifelse(sub_is_group_int, 1, ifelse(is_subgroup_int,1,sub_int.list.final[,3])),
                                                                          edge.color = ifelse(sub_is_group_int, ifelse(group_is_sub,'purple','black'), ifelse(is_subgroup_int,'green','grey')),
                                                                          knots = sub_int.list.final[,4],
                                                                          knot.size = ifelse(sub_is_group_int,3,ifelse(is_subgroup_int,3,1)),
                                                                          knot.color = ifelse(sub_is_group_int, ifelse(group_is_sub_int,'purple','black'), ifelse(is_subgroup_int,'green','grey')),
                                                                          knot.borders = TRUE,
                                                                          knot.border.width = 1,
                                                                          DoNotPlot=FALSE)
        } else {
          if(interaction_flag){
            output[['subgroup']][['interaction_fig']][[p]] = qgraph::qgraph(sub_int.list.final[,1:3],
                                                                            layout='circle',
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
      id_subpaths_all = id_subpaths + id_subpaths_all
       }
  }
  
  #### Create Group-Level Figure ####
  if(verbose){print('Creating plots for group.', quote = FALSE)}
  groupcut = output[["function_parameters"]][["groupcutoff"]] 
  temp_counts = output[["group"]][["group_paths_proportions"]] + id_subpaths_all # ID subgroup-level paths; zero if no subgroups
  if(output$function_parameters$subgroup && output$function_parameters$heuristic == 'GIMME'){
    temp_counts = temp_counts*ifelse(rowSums((moderated_group_by_subgroup < 0), dims=2),-1,1)
  }
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
  is_group = between(edge.list[,3], groupcut, 1)
  is_subgroup = edge.list[,3] > 1
  group_is_sub = sign(edge.list[,3]) < 0

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
                                                             edge.width = ifelse(is_group,1,ifelse(is_subgroup,1,edge.list[,3])),
                                                             edge.color = ifelse(is_group,ifelse(group_is_sub,'purple','black'),ifelse(is_subgroup,'green','grey')),
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
    is_group = between(int.list.final[,3], groupcut, 1)
    is_subgroup = int.list.final[,3] > 1
    group_is_sub = sign(int.list.final[,3]) < 0
    
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
                                                              edge.width = ifelse(is_group, 1, ifelse(is_subgroup,1,int.list.final[,3])),
                                                              edge.color = ifelse(is_group, ifelse(group_is_sub,'purple','black'), ifelse(is_subgroup,'green','grey')), 
                                                              knots = int.list.final[,4],
                                                              knot.size = ifelse(is_group,3,ifelse(is_subgroup,3,1)),
                                                              knot.color = ifelse(is_group, ifelse(group_is_sub,'purple','black'), ifelse(is_subgroup,'green','grey')),
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
                                                              edge.color = 'black',
                                                              knots = int.list.final[,4],
                                                              knot.size = 1,
                                                              knot.color = 'black',
                                                              knot.borders = TRUE,
                                                              knot.border.width = .5,
                                                              DoNotPlot=FALSE)
    }
    
  } else {
    if (interaction_flag) {
      output[['group']][['interaction_fig']] = qgraph::qgraph(int.list,
                                                              layout='circle',
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
