#' multiREG Model Selection
#' @param x Matrix of predictor variables.
#' @param y Vector of criterion variable.
#' @param selection_crit Model selection criteria (e.g., Bayesian information criterion, cross-validation).
#' @param alpha Elastic-net parameter.
#' @param penalty.factor User-set penalty matrix.
#' @return Returns best model parameters based on fit indices.
#' @keywords internal 
model_selection = function(x = NULL,
                           y = NULL,
                           selection_crit = NULL,
                           alpha = NULL,
                           penalty.factor = NULL){
  
  if (selection_crit == 'cv'){
    
    if (length(y) < 1000){ nfolds = round(length(y)/10) } else { nfolds = 100 }
    cvfit = glmnet::cv.glmnet(x = x,
                              y = y,
                              type.measure = 'mse',
                              nfolds = nfolds,
                              alpha = alpha,
                              penalty.factor = penalty.factor)
                    
    final_coefs = coef(cvfit, s = 'lambda.1se')
  } else {

    crit = selection_crit
    fit = glmnet::glmnet(x = x,
                         y = y,
                         alpha = alpha,
                         penalty.factor = penalty.factor)
    
    coef = coef(fit)
    lambda = fit$lambda
    df = fit$df

    yhat = cbind(1, x) %*% coef
    residuals = y - yhat
    mse = colMeans(as.array(residuals^2))
    sse = colSums(as.array(residuals^2))
    n = fit$nobs

    nvar = df + 1
    bic = n * log(mse) + nvar * log(n)
    aic = n * log(mse) + 2 * nvar
    aicc = aic + (2 * nvar * (nvar + 1)) / (n - nvar - 1)
    hqc = n * log(mse) + 2 * nvar * log(log(n))
    
    # tLL = fit$nulldev - glmnet::deviance(fit)
    # k = fit$df
    # n = fit$nobs
    # 
    # bic = log(n, base = exp(1)) * k - tLL
    # aic = 2 * k - tLL
    # aicc = (2 * k - tLL) + ((2 * k * k) + 2 * k)/(n - k - 1)
    # hqc = 2 * k * log(log(n, base = exp(1)), exp(1)) - tLL

    crit = switch(crit, bic=bic, aic=aic, aicc=aicc, hqc=hqc)
    selected = which(crit == min(crit))

    final_coefs = coef(fit, s = fit$lambda[selected])
       
  }
  return(final_coefs)
}
