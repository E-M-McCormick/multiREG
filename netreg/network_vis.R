#' @name network_vis
#' @aliases network_vis
#' @title Network Visualization following Regularized Search
#' @description This function utilizes output from network_reg to build network plots.
#' @usage 
#' network_vis(output = NULL)
#'
#'@param output Output object from network_reg
network_vis = function(output = NULL){
  library(qgraph)
  
  subnames = names(output[!names(output) %in% c('group','function_parameters','variablenames')])
  for (sub in subnames){
    temp_regmat = output[[sub]][['regression_matrix']]
    contemp = temp_regmat[!grepl('Lag', rownames(temp_regmat)) & !grepl('_x_', rownames(temp_regmat)),]
    lagged = temp_regmat[grepl('Lag', rownames(temp_regmat)) & !grepl('_x_', rownames(temp_regmat)),]
    interacted = temp_regmat[grepl('_x_', rownames(temp_regmat)),]
    
    con.list = cbind(which(contemp!=0 & !is.na(contemp), arr.ind=TRUE), contemp[contemp!=0 & !is.na(contemp)])
    lag.list = cbind(which(lagged!=0 & !is.na(lagged), arr.ind=TRUE), lagged[lagged!=0 & !is.na(lagged)])
    int.list = cbind(which(interacted!=0 & !is.na(interacted), arr.ind=TRUE), interacted[interacted!=0 & !is.na(interacted)])
    edge.list = rbind(con.list, lag.list)
    
    is_exogenous = rownames(contemp) %in% output[['variablenames']][['exogenous_vars']]
    is_lagged = grepl('Lag',rownames(edge.list)) & !grepl('_x_', rownames(edge.list))
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
      for (i in 1:nrow(int.list)){
        int.list2 = int.list1 = int.list
        rownames(int.list1) = sub('_x_.*','',rownames(int.list))
        int.list1[,1] = unique(con.list[rownames(con.list) %in% rownames(int.list1),1])
        rownames(int.list2) = sub('.*_x_','',rownames(int.list))
        int.list2[,1] = unique(con.list[rownames(con.list) %in% rownames(int.list2),1])
        int.list.final = rbind(int.list1,int.list2)
        
        is_lagged = grepl('Lag',rownames(int.list.final)) & !grepl('_x_', rownames(int.list.final))
        rownames(int.list.final)[is_lagged] = sub('Lag','',rownames(int.list.final)[is_lagged])
        is_pos = int.list.final[,3] > 0
        
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
      }
    }
  }
}
  
  
qgraph(edge.list,
       layout='circle',
       lty=ifelse(is_lagged,2,1),
       edge.labels = moderators,
       #edge.color = colors,
       posCol='red',
       negCol='blue',
       esize = 5,
       edge.label.position = ifelse(is_interact,label_pos,0),
       edge.label.bg = F,
       edge.label.cex = .5,
       curveAll = F,
       curveDefault = 1,
       cut = .75,
       parallelEdge = T,
       fade = F,
       labels = sub('_',' ', rownames(contemp)),
       label.cex = 1,
       shape = ifelse(is_exogenous,'square','circle'),
       DoNotPlot=F)
  
  colors = ifelse(edge.list[,3] > 7,'black','grey')
  pdf(file.path(getwd(), "summaryPathsPlot.pdf"))
  plot(fig)
  dev.off()

  edge.list = rbind(con.list, lag.list, int.list)
  
  is_interact = grepl('_x_', rownames(edge.list))
  is_exogenous = rownames(contemp) %in% output[['variablenames']][['exogenous_vars']]
  
  moderators = ifelse(is_interact,sub('.*_x_','',rownames(edge.list)),'')
  rownames(edge.list)[is_interact] = sub('_x_.*','',rownames(edge.list)[is_interact])
  edge.list[,1][is_interact] = edge.list[!is_interact,1][names(edge.list[,1][is_interact])]
  is_lagged = grepl('Lag',rownames(edge.list)) & !grepl('_x_', rownames(edge.list))
  rownames(edge.list)[is_lagged] = sub('Lag','',rownames(edge.list)[is_lagged])
  label_pos = 0.5+(runif(length(edge.list))-0.5)/2
  output[[sub]][['main_effects_fig']] = fig = qgraph(edge.list,
                                                     layout='circle',
                                                     lty=ifelse(is_lagged,2,1),
                                                     posCol='red',
                                                     negCol='blue',
                                                     esize = 5,
                                                     edge.labels = moderators,
                                                     edge.label.position = ifelse(is_interact,label_pos,0),
                                                     edge.label.bg = F,
                                                     edge.label.cex = .5,
                                                     parallelEdge = T,
                                                     fade = F,
                                                     labels = sub('_',' ', rownames(contemp)),
                                                     label.cex = 1,
                                                     shape = ifelse(is_exogenous,'square','circle'),
                                                     DoNotPlot=F)