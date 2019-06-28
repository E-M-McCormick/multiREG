#' @name manage_output
#' @title Manage Output Parameters Following Regularization with LASSO
#' @keywords internal 
manage_output = function(out = NULL, plot = NULL, output = NULL){
  if (!is.null(out)){
    yvarnames = output[['variablenames']][['y_vars']]
    groupcutoff = output[['function_parameters']][['groupcutoff']]
    
    # Create Output Directory if Needed
    if (!dir.exists(out)){
      print('Creating output directories', quote = FALSE)
      dir.create(out);
    }
    print('Writing output to file.', quote = FALSE)
    
    # Write out Function Arguments & Variable Names
    capture.output(print('Function Arguments', quote = FALSE), 
                   print(output$function_parameters), 
                   print('Variable Names', quote = FALSE), 
                   print(output$variablenames), 
                   file = paste(out, 'function_summary.txt', sep=.Platform$file.sep))
    
    # Write Group Level Data with Plots if Needed
    write.csv(output$group$group_paths_counts[, colnames(output$group$group_paths_counts) %in% yvarnames],
              file = paste(out, 'groupPathCountsMatrix.csv', sep=.Platform$file.sep))
    write.csv(output$group$group_paths_present[, colnames(output$group$group_paths_present) %in% yvarnames],
              file = paste(out, 'groupPathCountsPresent.csv', sep=.Platform$file.sep))
    if (plot) {
      pdf(file.path(out, 'groupMainEffectsPlots.pdf'))
      plot(output[['group']][['main_effects_fig']])
      dev.off()
      pdf(file.path(out, 'groupInteractionsPlots.pdf'))
      plot(output[['group']][['interaction_fig']])
      dev.off()
    }
    
    # Write Individual Level Data with Plots if Needed
    dir.create(paste(out, 'individual', sep=.Platform$file.sep))
    indpaths = data.frame()
    pathtypes = output[['group']][['group_paths_proportions']]
    pathtypes[pathtypes >= groupcutoff] = 'group'
    pathtypes[!is.na(pathtypes) & pathtypes > 0 & pathtypes != 'group'] = 'individual'
    pathtypes[pathtypes != 'group' & pathtypes != 'individual'] = 'none'
    
    subnames = names(output[!names(output) %in% c('group','function_parameters','variablenames')])
    for (sub in subnames){
      write.csv(output[[sub]][['regression_matrix']][, colnames(output$group$group_paths_counts) %in% yvarnames],
                file = paste(out, 'individual', paste0(sub,'Betas.csv'), sep=.Platform$file.sep))
      if (plot) {
        pdf(file.path(out, 'individual', paste0(sub,'MainEffectsPlot.pdf')))
        plot(output[[sub]][['main_effects_fig']])
        dev.off()
        pdf(file.path(out, 'individual', paste0(sub,'InteractionsPlot.pdf')))
        plot(output[[sub]][['interaction_fig']])
        dev.off()
      }
      ind=cbind(which(output[[sub]][['regression_matrix']] != 0 & !is.na(output[[sub]][['regression_matrix']]), arr.ind = TRUE, useNames = F), 
                output[[sub]][['regression_matrix']][output[[sub]][['regression_matrix']] != 0 & !is.na(output[[sub]][['regression_matrix']])])
      temp = ind
      temp[,1]=rownames(output$group$group_paths_counts)[ind[,1]]
      temp[,2]=rownames(output$group$group_paths_counts)[ind[,2]]
      temp = cbind(sub, temp)
      temp = cbind(temp, rep(0, nrow(temp)))
      for (r in 1:nrow(temp)){ temp[r,5] = pathtypes[temp[r,2], temp[r,3]] }
      colnames(temp) = c('file','iv','dv','beta_estimate','type')
      temp = temp[order(temp[,'type']),]
      indpaths = rbind(indpaths, temp)
    }
    write.csv(indpaths,
              file = paste(out, 'indivPathEstimates.csv', sep=.Platform$file.sep))
    setwd(out)
  }
}
