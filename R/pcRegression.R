#' pcRegression
#'
#' @description pcRegression does a linear model fit of principal components
#' and a batch (categorical) variable
#' @param pca.data a list as created by 'prcomp', pcRegression needs
#' \itemize{
#'  \item \code{x} -
#' the principal components (PCs, correctly: the rotated data) and
#' \item \code{sdev} - the standard deviations of the PCs)
#' }
#' @param batch vector with the batch covariate (for each cell)
#' @param n_top the number of PCs to consider at maximum
#' @param tol truncation threshold for significance level, default: 1e-16
#' @return List summarising principal component regression
#' \itemize{
#'  \item \code{maxVar} - the variance explained by principal component(s)
#'  that correlate(s) most with the batch effect
#'  \item \code{PmaxVar} - p-value (returned by linear model) for the
#'  respective principal components (related to \code{maxVar})
#'  \item \code{pcNfrac} - fraction of significant PCs among the \code{n_top} PCs
#'  \item \code{pcRegscale} - 'scaled PC regression', i.e. total variance of PCs which correlate significantly with batch covariate (FDR<0.05) scaled by the total variance of \code{n_top} PCs
#'  \item \code{maxCorr} - maximal correlation of \code{n_top} PCs with batch covariate
#'  \item \code{maxR2} - maximal coefficient of determination of \code{n_top} PCs with batch covariate
#'  \item \code{msigPC} - scaled index of the smallest PC that correlates significantly with batch covariate (FDR<0.05), i.e. \code{msigPC=1} if PC_1 is significantly correlated with the batch covariate and \code{msigPC=0} if none of the \code{n_top} PCs is significantly correlated
#'  \item \code{maxsigPC} - similar to \code{msigPC}, scaled index of the PC with maximal correlation of \code{n_top} PCs with batch covariate
#'  \item \code{R2Var} - sum over Var(PC_i)*r2(PC_i and batch) for all i
#' \item \code{ExplainedVar} - explained variance for each PC
#' \item \code{r2} - detailed results of correlation (R-Square) analysis
#' }
#' @examples
#'     testdata <- create_testset_multibatch(n.genes=1000, n.batch=3, plattform='any')
#'     pca.data <- prcomp(testdata$data, center=TRUE)
#'     pc.reg.result <- pcRegression(pca.data, testdata$batch)
#' @importFrom stats t.test lm aov p.adjust
#' @export
pcRegression <- function(pca.data, batch,n_top=50, tol=1e-16){
  batch.levels <- unique(batch)
#make sure you do not try to assess more PCs than actually computed
  pca_rank = ncol(pca.data$x)
  max_comps <- min(pca_rank, n_top)
  if (length(pca.data$sdev) > pca_rank) {
    pca.data$sdev = pca.data$sdev[1:pca_rank]
    }

  if (length(batch.levels) == 2) {
    #r2.batch.raw <- r2.batch

    # for-loop replaced by correlate.fun and apply
    r2.batch <- apply(pca.data$x, 2, correlate.fun_two, batch, batch.levels)
    r2.batch <- t(r2.batch)
    colnames(r2.batch) <- c('R.squared', 'p.value.lm', 'p.value.t.test')

    r2.batch[r2.batch[,2] < tol, 2] <- tol
    r2.batch[r2.batch[,3] < tol, 3] <- tol
  } else {


  #r2.batch.raw <- r2.batch

  r2.batch <- apply(pca.data$x, 2, correlate.fun_gen, batch)

  r2.batch <- t(r2.batch)
  colnames(r2.batch) <- c('R.squared', 'p.value.lm', 'p.value.F.test')
  # for-loop replaced by correlate.fun and apply
  #for (k in 1:dim(r2.batch)[1]){
  #    a <- lm(pca.data$x[,k] ~ batch)
  #    r2.batch[k,1] <- summary(a)$r.squared #coefficient of determination
  #    r2.batch[k,2] <- summary(a)$coefficients['batch',4] #p-value (significance level)
  #}

  r2.batch[r2.batch[,2] < tol, 2] <- tol
  r2.batch[r2.batch[,3] < tol, 3] <- tol
  }

  argmin <- which(r2.batch[, 2] == min(r2.batch[, 2]))
  normal <- sum(pca.data$sdev^2)
  var <- round((pca.data$sdev)^2 / normal * 100,1)
  batch.var <- sum(r2.batch[,1]*var)/100
  setsignif <- p.adjust(r2.batch[1:max_comps,2], method = 'BH') < 0.05
  pcCorr <- sqrt(r2.batch[1:max_comps, 1])
  result <- list()
  result$maxVar <- var[argmin]
  result$PmaxVar <- r2.batch[argmin,2]
  result$pcNfrac <- mean(setsignif)
  result$pcRegscale <- sum(var[1:max_comps][setsignif])/sum(var[1:max_comps])
  result$maxCorr <- max(pcCorr)
  result$maxR2 <- max(r2.batch[1:max_comps, 1])
  result$msigPC   <- 1 - (min(c(which(setsignif), max_comps + 1)) - 1) / max_comps
  result$maxsigPC <- 1 - (min(c(which(pcCorr == max(pcCorr[setsignif])), max_comps + 1)) - 1) / max_comps
  result$R2Var <- batch.var
  result$ExplainedVar <- var
  result$r2 <- r2.batch
  result
}


