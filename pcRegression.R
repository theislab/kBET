#this function does a linear model fit of principal components and a batch (categorical) variable
# pca.data: a list as created by 'prcomp', pcRegression needs
#   * $x: the principal components (PCs, correctly: the rotated data)
#   * $sdev: the standard deviations of the PCs
# batch: vector with the batch covariate (for each cell)
# tol: truncation threshold for significance level, default: 1e-16

pcRegression <- function(pca.data, batch, tol=1e-16){
  batch.levels <- unique(batch) 
  if(length(batch.levels)==2){
    r2.batch <- matrix(NA,nrow=dim(pca.data$x)[2], ncol=3)
    colnames(r2.batch) <- c('R.squared', 'p.value.lm', 'p.value.t.test')
    #r2.batch.raw <- r2.batch
    for (k in 1:dim(r2.batch)[1]){
      a <- lm(pca.data$x[,k] ~ batch)
      r2.batch[k,1] <- summary(a)$r.squared
      r2.batch[k,2] <- summary(a)$coefficients['batch',4]
      t.test.result <- t.test(pca.data$x[batch==batch.levels[1],k], 
                              pca.data$x[batch==batch.levels[2],k], paired = FALSE)
      r2.batch[k,3] <- t.test.result$p.value
      
      #a <- lm(testset[,i] ~ batch)
      # r2.batch.raw[i] <- summary(a)$r.squared
    }
    r2.batch[r2.batch[,2]<tol,2] <- tol
    r2.batch[r2.batch[,3]<tol,3] <- tol
  }else{
  r2.batch <- matrix(NA,nrow=dim(pca.data$x)[2], ncol=2)
  colnames(r2.batch) <- c('R.squared', 'p.value.lm')
  #r2.batch.raw <- r2.batch
  correlate.fun <- function(rot.data, batch){
      a <- lm(rot.data[,k] ~ batch)
      result <- numeric(2)
      result[1] <- summary(a)$r.squared #coefficient of determination
      result[2] <- summary(a)$coefficients['batch',4] #p-value (significance level)
  }
  r2.batch <- apply(pca.data$x, 2, correlate.fun, batch)
  # for-loop replaced by correlate.fun and apply
  #for (k in 1:dim(r2.batch)[1]){
  #    a <- lm(pca.data$x[,k] ~ batch)
  #    r2.batch[k,1] <- summary(a)$r.squared #coefficient of determination
  #    r2.batch[k,2] <- summary(a)$coefficients['batch',4] #p-value (significance level)
  #}
  
  r2.batch[r2.batch[,2]<tol,2] <- tol
  
  }
  
  argmin <- which(r2.batch[,2]==min(r2.batch[,2]))
  normal <- sum(pca.data$sdev^2)
  var <- round((pca.data$sdev)^2 / normal *100,1)
  batch.var <- sum(r2.batch[,1]*var)/100
  result <- list()
  result$maxVar <- var[argmin]
  result$PmaxVar <- r2.batch[argmin,2]
  result$R2Var <- batch.var
  result$ExplainedVar <- var
  result$r2 <- r2.batch
  return(result)
}
