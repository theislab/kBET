#this function does a linear model fit of principal components and a batch (categorical) variable
# pca.data: a list as created by 'prcomp', pcRegression needs
#   * $x: the principal components (PCs, correctly: the rotated data)
#   * $sdev: the standard deviations of the PCs
# batch: vector with the batch covariate (for each cell)
# tol: truncation threshold for significance level, default: 1e-16

pcRegression <- function(pca.data, batch, tol=1e-16){
  batch.levels <- unique(batch)
  
  
  if(length(batch.levels)==2){
    #r2.batch.raw <- r2.batch
    correlate.fun <- function(rot.data, batch){
        a <- lm(rot.data ~ batch)
        result <- numeric(2)
        result[1] <- summary(a)$r.squared #coefficient of determination
        result[2] <- summary(a)$coefficients['batch',4] #p-value (significance level)
        t.test.result <- t.test(rot.data[batch==batch.levels[1]],
                                rot.data[batch==batch.levels[2]], paired = FALSE)
        result[3] <- t.test.result$p.value
        return(result)
    }
    # for-loop replaced by correlate.fun and apply
    r2.batch <- apply(pca.data$x, 2, correlate.fun, batch)
    r2.batch <- t(r2.batch)
    colnames(r2.batch) <- c('R.squared', 'p.value.lm', 'p.value.t.test')
    
    r2.batch[r2.batch[,2]<tol,2] <- tol
    r2.batch[r2.batch[,3]<tol,3] <- tol
  }else{

 
  #r2.batch.raw <- r2.batch
  correlate.fun <- function(rot.data, batch){
      a <- lm(rot.data ~ batch)
      result <- numeric(2)
      result[1] <- summary(a)$r.squared #coefficient of determination
      result[2] <- summary(a)$coefficients['batch',4] #p-value (significance level)
      return(result)
  }
  r2.batch <- apply(pca.data$x, 2, correlate.fun, batch)
  
  r2.batch <- t(r2.batch)
  colnames(r2.batch) <- c('R.squared', 'p.value.lm')
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
