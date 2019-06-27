# multiLASSO Model Selection

dev = deviance(fit)
bic= log(900, base = exp(1))*fit$df + deviance(fit)
plot(bic)


x = subset_predictors
y = tempdata[, colnames(tempdata) %in% varname]
#crit=match.arg(crit)
crit = 'bic'
n=length(y)
model = glmnet(x = x, 
               y = y, 
               alpha = alpha,
               penalty.factor = initial_penalties[!colnames(tempdata) %in% varname, varname])













coef = coef(model)
lambda = model$lambda
df = model$df
yhat=cbind(1,x)%*%coef
residuals = (y - yhat)
mse = colMeans(residuals^2)
sse = colSums(residuals^2)
nvar = df + 1
bic = n*log(mse)+nvar*log(n)
aic = n*log(mse)+2*nvar
aicc = aic+(2*nvar*(nvar+1))/(n-nvar-1)
hqc = n*log(mse)+2*nvar*log(log(n))
sst = (n-1)*var(y)
r2 = 1 - (sse/sst)
adjr2 = (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1))
