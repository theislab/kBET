
#' create_testset_multibatch
#'
#' @description \code{create_testset_multibatch} simulates single-cell RNA-seq data
#' as zero-inflated negative binomial counts. The observe several batches that differ from each other - the difference is
#' randomly chosen from the Beta distribution with a = 1, b = 9, i.e. the expected fraction of varied genes is 1/10.
#' @param n.genes nummber of sampled genes
#' @param n.batch number of batches in the data set
#' @param plattform maximum number of cells (samples) per batch. There are three modes:
#'  \enumerate{
#'    \item \code{C1} - 96 cells (max)
#'    \item \code{dropseq} - 1000 cells (max)
#'    \item \code{any} - 300 cells (max)
#'    }
#'  The number of cell per batch is varied by a Beta distributed random variable to simulate a quality control, where
#'  80\% of the samples pass.
#'
#' @return list object
#' \enumerate{
#'    \item \code{data} - simulated data set with #rows: cells, columns: 'genes'
#'    \item \code{batch} - vector with batch identity for each cell
#'    \item \code{comment} - 'testSet_multibatch_<n.genes>genes_<n.batch>batches_<plattform>'
#'    }
#'
#' @importFrom MASS rnegbin
#' @importFrom stats rbeta rgamma rbinom
#' @export
#'
#' @name create_testset_multibatch
#' @examples
#'     testdata <- create_testset_multibatch(n.genes=1000, n.batch=3, plattform='any')
#'
create_testset_multibatch <- function(n.genes=1000, n.batch=3, plattform = c('C1', 'dropseq', 'any'))
{

  #let us simulate several plattforms:
  #C1 - 96 cells (max)
  #dropseq - 3000 cells (max)
  #any - 500 cells (max)
  #we vary the amount of cells per batch to simulate the variance from the QC

  #preliminary initialisations:
  alpha <- 0.05
  mu<- rbeta(n.genes-1, 2,5)*100
  mu <- c(mean(mu), mu)
  b0 <- -1.5

  cells <- data.frame(plattform = c('C1', 'dropseq', 'any'), batch.size = c(96, 1500, 300))
  #get sample size per batch
  samples <- rbeta(n.batch,8,2)*cells$batch.size[cells$plattform==plattform]

  #simulate means of the all batches
  #we choose the variability between batches with the beta distribution. Expected value= a/(a+b)= 1/10 here.
  mu.batch <- sapply(floor(n.genes* rbeta(n.batch, 1,9)),
                     function(x, mu){c(rgamma(x,1), rep(1, n.genes-x))*mu}, mu)
  #model drop-out rate as logistic model
  b2 <- apply(mu.batch, 2, function(x){1/quantile(x,0.5)})
  decay.prob2 <- apply(mu.batch,2, function(mu, b0){
    b1 <- 1/quantile(mu,0.5)
    res <- 1/(1+exp(-(b0+b1*mu)))}, b0)

  #simulate counts in all batches
  testset <- sapply(seq_len(n.genes),function(k,sample.size, decay.prob, mu){
    unlist(sapply(seq_len(n.batch), function(x,sample.size, decay.prob, mu){rnegbin(sample.size[x], mu=mu[k,x], 1)*
        rbinom(sample.size[x], 1, decay.prob[k,x])},
        sample.size, decay.prob, mu))}, samples, decay.prob2, mu.batch)

  result <- list()

  result$data <- testset #rows: cells, columns: 'genes'
  result$batch <- unlist(sapply(1:n.batch, function(x,y){rep(x, y[x])}, samples))
  result$comment <- paste0('testSet_multibatch_', n.genes, 'genes_', n.batch, 'batches_', plattform)

 return(result)
}
