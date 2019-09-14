library(mvtnorm) 
library(MASS)
library(urca)

set.seed(1107)

## Simulation Condition for factors = No. of Variable and % of Sparsity
#condition.matrix <- matrix(c(20, 20, 200, 200, .05, .1, .05, .1), 4, 2)
#colnames(condition.matrix)  <- c("No.Var", "Sparsity")
condition.matrix <- matrix(c(20, 200), 2, 1)
colnames(condition.matrix)  <- c("No.Var")
## Simulation Condition for factor No. of Time Points (will cross within condition.matrix)
No.TP <- c(60, 200, 1000) # Separate No. of TP conditions from No. of VAR and % of Sparsity to avoid resampling and keeping path assignment consistent
## Simulation Condition for factor Path Strength (not full cross, mixed within condition.matrix, with fix proportion 1:4:5)
Strength.S <- c(.6, .3, .1) # strong, moderate, weak path weight in small network
Strength.L <- c(.3, .1, .05) # for the large and sparse network
#Strength.LD <- c(.1, .05, .01) # for the large dense network
# Strength <- c(.1, .1, .05) # even 50/50 .1 and .05 does NOT worked for the large dense network
# Strength <- c(.1, .05, .01) 
## Data generating function (12 .Rdata will be generated, called "SimDataCond_XX.YY.ZZ.RData")
genData <- function(condition.matrix = NULL, n = 1000){ 
  
  datatemp <- list() # store individual dataset of each condition into a list (4 in total)
  
  for(i in 2:nrow(condition.matrix)){  
    ## Create Parameter Matricies
    A <- matrix(0, nrow=condition.matrix[i,"No.Var"], ncol=condition.matrix[i,"No.Var"])
    Phi <- matrix(0, nrow=condition.matrix[i,"No.Var"], ncol=condition.matrix[i,"No.Var"])
    Psi <- diag(condition.matrix[i,"No.Var"])
    ## Create a vector to store diagonal index (no diagonal assignment for A matrix)
    x<-c(1:ncol(A))
    diag.list <- ncol(A)*(x-1)+x
    # Alist <- c(1:400)[-c(diag.list)]
    size = ncol(A)*nrow(A)*.05
    ## Random sampling path locations
    Pathposition <- sample(c(1:(ncol(A)*nrow(A))), size=size)
    Apath1 <- vector()
    Apath2 <- vector()
    Apath3 <- vector()
    Apathall <- vector()
    Phipath1 <- vector()
    Phipath2 <- vector()
    Phipath3 <- vector()
    Phipathall <- vector()
    pathlist <- list()
    if(condition.matrix[i,"No.Var"]==20){ 
      ## Random order of path stength (3 levels) with a proportion of strong:moderate:weak = 1:4:5
      Weightorder <- sample(c(Strength.S[1], Strength.S[2], Strength.S[3]), size=size, prob=c(1,4,5), replace = TRUE)
      for(r in 1:size){
        # if the random position happens to be on the diagonal, it will be assigned to Phi (as AR)
        if(Pathposition[r] %in% diag.list){Phi[Pathposition[r]] <- Weightorder[r]}
        # if not, it has 50/50 chance of being assigned to A or Phi
        else{
          A.Phi.Assign <- sample(c(0,1), size=1)
          if(A.Phi.Assign==0) {A[Pathposition[r]]=Weightorder[r]} 
          if(A.Phi.Assign==1) {Phi[Pathposition[r]]=Weightorder[r]}
        }
      }
      # store the position of paths in each matrix
      Apath1<-which(A==Strength.S[1])
      Apath2<-which(A==Strength.S[2])
      Apath3<-which(A==Strength.S[3])
      Apathall <- which(A!=0)
      Phipath1<-which(Phi==Strength.S[1])
      Phipath2<-which(Phi==Strength.S[2])
      Phipath3<-which(Phi==Strength.S[3])
      Phipathall <- which(Phi!=0)
      pathlist <- list(Apath1, Apath2, Apath3, Apathall,
                       Phipath1, Phipath2, Phipath3, Phipathall)
      save(pathlist,file="PathLoc_V20.RData")
    }
    else{
      Weightorder <- sample(c(Strength.L[1], Strength.L[2], Strength.L[3]), size=size, prob=c(1,4,5), replace = TRUE)
      for(r in 1:size){
        # if the random position happens to be on the diagonal, it will be assigned to Phi (as AR)
        if(Pathposition[r] %in% diag.list){Phi[Pathposition[r]] <- Weightorder[r]}
        # if not, it has 50/50 chance of being assigned to A or Phi
        else{
          A.Phi.Assign <- sample(c(0,1), size=1)
          if(A.Phi.Assign==0) {A[Pathposition[r]]=Weightorder[r]} 
          if(A.Phi.Assign==1) {Phi[Pathposition[r]]=Weightorder[r]}
        }
      }
      # store the position of paths in each matrix
      Apath1<-which(A==Strength.L[1])
      Apath2<-which(A==Strength.L[2])
      Apath3<-which(A==Strength.L[3])
      Apathall <- which(A!=0)
      Phipath1<-which(Phi==Strength.L[1])
      Phipath2<-which(Phi==Strength.L[2])
      Phipath3<-which(Phi==Strength.L[3])
      Phipathall <- which(Phi!=0)
      pathlist <- list(Apath1, Apath2, Apath3, Apathall,
                       Phipath1, Phipath2, Phipath3, Phipathall)
      save(pathlist,file="PathLoc_V200.RData")
    }
    # generate individual time series data (No.Var*No.TP)
    for(t in 1:3){
      for (p in 1:n){
        negA1 <- solve(diag(condition.matrix[i,"No.Var"])-A)
        time <- matrix(0, nrow = condition.matrix[i,"No.Var"], ncol = No.TP[t] +2*No.TP[t])
        time1 <- matrix(0, nrow = condition.matrix[i,"No.Var"], ncol = No.TP[t]+2*No.TP[t])
        noise <- negA1 %*% t(mvrnorm(n = condition.matrix[i,"No.Var"]*(No.TP[t]+2*No.TP[t]),
                                     rep(0,condition.matrix[i,"No.Var"]), Psi, empirical = TRUE))
        time[,1] <- noise[,1]
        time1[,1] <- negA1 %*% Phi %*% time[,1] + noise[,1]
        time[,2] <- time1[,1]
        for (j in 2:(No.TP[t]+2*No.TP[t])){
          time1[,j]  <- negA1 %*% Phi %*% time[,(j)] + noise[,j]
          if (j<(No.TP[t]+2*No.TP[t]))
            time[,(j+1)] <- time1[,j]
        }
        datatemp[[p]] <- t(time[,(2*No.TP[t]):(3*No.TP[t])])
        colnames(datatemp[[p]]) <- c(1:condition.matrix[i,"No.Var"])
        names(datatemp)[p]<- paste0('ind', p)
        V10 <- sample(c(1:condition.matrix[i,"No.Var"]), size=10)
        jotest <- ca.jo(datatemp[[p]][,V10], type="eigen", K=2, ecdet="none", spec="longrun")
        ## Weakly stantionary test by Johansen'test for VAR
        while(jotest@teststat[1]<jotest@cval[1,2]){
          time <- matrix(0, nrow = condition.matrix[i,"No.Var"], ncol = No.TP[t] +2*No.TP[t])
          time1 <- matrix(0, nrow = condition.matrix[i,"No.Var"], ncol = No.TP[t]+2*No.TP[t])
          noise <- negA1 %*% t(mvrnorm(n = condition.matrix[i,"No.Var"]*(No.TP[t]+2*No.TP[t]),
                                       rep(0,condition.matrix[i,"No.Var"]), Psi, empirical = TRUE))
          time[,1] <- noise[,1]
          time1[,1] <- negA1 %*% Phi %*% time[,1] + noise[,1]
          time[,2] <- time1[,1]
          for (j in 2:(No.TP[t]+2*No.TP[t])){
            time1[,j]  <- negA1 %*% Phi %*% time[,(j)] + noise[,j]
            if (j<(No.TP[t]+2*No.TP[t]))
              time[,(j+1)] <- time1[,j]
          }
          datatemp[[p]] <- t(time[,(2*No.TP[t]):(3*No.TP[t])])
          colnames(datatemp[[p]]) <- c(1:condition.matrix[i,"No.Var"])
          names(datatemp)[p]<- paste0('ind', p)
          V10 <- sample(c(1:condition.matrix[i,"No.Var"]), size=10)
          jotest <- ca.jo(datatemp[[p]][,V10], type="eigen", K=2, ecdet="none", spec="longrun")
        }
        condition.factor <- levels(interaction(condition.matrix[i,"No.Var"], No.TP[t]))
        save(datatemp,file=paste("SimDataCond_", condition.factor, ".RData", sep = ""))
      }
    }
  }
}

genData(condition.matrix = condition.matrix, n=1000)





