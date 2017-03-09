#' kBET - k-nearest neighbour batch effect test
#'
#' @description kBET runs a chi square test to evaluate the probability of a batch effect.
#'
#' @param df dataset (rows: samples, columns: parameters)
#' @param batch batch id for each cell or a data frame with both condition and replicates
#' @param k0 number of nearest neighbours to test on (neighbourhood size)
#' @param knn a set of nearest neighbours for each cell (optional)
#' @param testSize number of data points to test, (10 percent sample size default)
#' @param heuristic compute an optimal neighbourhood size k
#' @param stats to create a statistics on batch estimates, evaluate 'stats' subsets
#' @param alpha significance level
#' @param addTest perform an LRT-approximation to the multinomial test AND a multinomial exact test (if appropriate)
#' @param plot if stats > 10, then a boxplot of the resulting rejection rates is created
#' @return List with test summary - a rejection rate for the data, an expected rejection rate for random labeling and the significance for the observed result
#' @examples
#' batch.estimate <- kBET(data,batch)

#' @importFrom FNN get.knn
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @include bisect.R
#' @include kBET-utils.R
#' @name kBET
#' @export
kBET <- function(df, batch, k0=NULL,knn=NULL, testSize=NULL,heuristic=FALSE, stats=100, alpha=0.05, addTest = FALSE, verbose=TRUE, plot = TRUE){


  if (plot==TRUE & heuristic==TRUE){
    #source('addalpha.R')
    colorset <- addalpha(brewer.pal(8, 'Set2'), 0.5)
    colorset.nt <- brewer.pal(8, 'Set2')
  }

  #create a subsetting mode:
  #if (is.data.frame(batch) && dim(batch)[2]>1){
  #  cat('Complex study design detected.\n')
  #  design.names <- colnames(batch)
  #  cat(paste0('Subset by ', design.names[[2]], '\n'))
  #  bio <- unlist(unique(batch[[design.names[2]]]))
    #drop subset batch
  #  new.batch <- base::subset(batch, batch[[2]]== bio[1])[,!design.names %in% design.names[[2]]]
  #  sapply(df[batch[[2]]==bio[1],], kBET, new.batch, k0, knn, testSize, heuristic,stats, alpha, addTest, verbose)

  #}


  #preliminaries:
  dof <- length(unique(batch))-1 #degrees of freedom
  if(is.factor(batch)){
    batch <- droplevels(batch)
  }
  frequencies <- table(batch)/length(batch)
  batch.shuff <- batch[sample.int(length(batch))]

  class.frequency <- data.frame(class = names(frequencies),
                                freq = as.numeric(frequencies))
  dataset <- df
  dim.dataset <- dim(dataset)
  stopifnot(class(stats) == 'numeric', stats>0)

  if (is.null(k0) || k0>=dim.dataset[1]){
    if(!heuristic){
      #default environment size: quarter the size of the largest batch
      k0=floor(mean(class.frequency$freq)*dim.dataset[1]/4)
    }else{
      #default environment size: size of the largest batch
      k0=floor(mean(class.frequency$freq)*dim.dataset[1])
    }

    cat('Initial neighbourhood size is set to ')
    cat(paste0(k0, '.\n'))
  }
  # find KNNs
  if (is.null(knn)){
    if (verbose) {
      cat('finding knns...')
      tic <- proc.time()
    }
    knn <- get.knn(dataset, k=k0, algorithm = 'cover_tree')
    if (verbose) {
      cat('done. Time:\n')
      print(proc.time() - tic)
    }
  }
  #set number of tests
  if (is.null(testSize) || dim.dataset[1]< testSize){
    test.frac <- 0.1
    testSize <- ceiling(dim.dataset[1]*test.frac)
    if (verbose) {
      cat('Number of kBET tests is set to ')
      cat(paste0(testSize, '.\n'))
    }
  }

  if(heuristic==TRUE){
    #source('bisect.R')
    myfun <- function(x,df,batch, knn){
      res <- kBET(df=df, batch=batch, k0=x, knn=knn, testSize=NULL, heuristic=FALSE, stats=10, alpha=0.05, addTest = FALSE, plot=FALSE, verbose=FALSE)
      result <- res$summary
      result <- result$kBET.observed[1]
    }
    #btw, when we bisect here that returns some interval with the optimal neihbourhood size
    if (verbose){
      cat('Determining optimal neighbourhood size ...')
    }
    opt.k <- bisect(myfun, bounds=c(10,k0), known=NULL, df, batch, knn)
    #result
    if(length(opt.k)>1){
      k0 <- opt.k[2]
      if(verbose==TRUE){
        cat('done.\n')
        #cat('The optimal neighbourhood size is determined.\n')
        cat('Size of neighbourhood is set to ')
        cat(paste0(k0, '.\n'))
      }
    }else{
      if(verbose==TRUE){
        cat('done.\n')
        #cat('The optimal neighbourhood size is determined.\n')
        cat('Heuristic did not change the neighbourhood.\n If results appear inconclusive, increase k0.')
        cat(paste0(k0, '.\n'))
      }
    }


    #heuristic was updated to an optimization based method instead of scanning method
    #plot of different test values now obsolete

    # if (plot==TRUE){
    #   k.frac <- 100*k/dim(df)[1]
    #   data.null <- matrix(unlist(heuristics['kBET.expected',]), ncol=test.k)
    #   data.kdep <- data.frame(neighbourhood=rep(c(k.frac, rev(k.frac)), 2),
    #                           result=c(data.heu[2,], rev(data.heu[4,]), #chi2 observed
    #                                    data.null[2,], rev(data.null[4,])), #chi2 expected
    #                           test= as.factor(rep(c('kBET','chi2 expected'), each=2*length(k))))
    #   data.kdep.points <- data.frame(neighbourhood=rep(k.frac,2), result= c(data.heu[1,], data.null[1,]),
    #                                  test= as.factor(rep(c('kBET','chi2 expected'), each=length(k))))
    #
    #   p.heu <- ggplot(data.kdep, aes(x=neighbourhood, y=result)) +
    #     geom_polygon(data=data.kdep, mapping=aes(x=neighbourhood, y=result, group=test, fill=test)) +
    #     geom_line(data=data.kdep.points, mapping=aes(x=neighbourhood, y=result, color=test)) +
    #     scale_color_manual(values = c('kBET' = colorset[1], 'chi2 expected'= colorset[2]),
    #                        name = 'Test', breaks=c('kBET','chi2 expected'),
    #                        labels = c('kBET', 'random\nassignment')) +
    #     scale_fill_manual(values = c('kBET' = colorset[1], 'chi2 expected'= colorset[2]),
    #                       name = 'Test', breaks=c('kBET','chi2 expected'),
    #                       labels = c('kBET', 'random\nassignment')) +
    #     labs(x= 'Neighbourhood size (in % sample size)', y = 'Rejection rate')+
    #     geom_vline(aes(xintercept = k.frac[opt.k]),
    #                linetype="dotted" , size=0.5)+
    #     theme_bw()
    #   print(p.heu)
    # }

  }

  if(addTest==TRUE){
    #initialize result list
  rejection <- list()
  rejection$summary <- data.frame(kBET.expected = numeric(4),
                                  kBET.observed = numeric(4),
                                  lrt.expected = numeric(4),
                                  lrt.observed = numeric(4)
  )
  kBET.expected <- numeric(stats)
  kBET.observed <- numeric(stats)
  kBET.signif <- numeric(stats)
  rejection$results   <- data.frame(tested = numeric(dim.dataset[1]),
                                    kBET.pvalue.test = rep(0,dim.dataset[1]),
                                    kBET.pvalue.null = rep(0, dim.dataset[1]),
                                    lrt.pvalue.test = rep(0,dim.dataset[1]),
                                    lrt.pvalue.null = rep(0, dim.dataset[1]))
  lrt.expected <- numeric(stats)
  lrt.observed <- numeric(stats)
  lrt.signif <- numeric(stats)
  #decide to perform exact test or not
  if (choose(k0+dof,dof)<5e5 & k0 <= min(table(batch))){
    exact.expected <- numeric(stats)
    exact.observed <- numeric(stats)
    exact.signif <- numeric(stats)

    rejection$summary$exact.expected <- numeric(4)
    rejection$summary$exact.observed <- numeric(4)
    rejection$results$exact.pvalue.test <- rep(0, dim.dataset[1])
    rejection$results$exact.pvalue.null <- rep(0, dim.dataset[1])
  }

  for (i in 1:stats){
    # choose a random sample from dataset (rows: samples, columns: parameters)
    idx.runs <- sample.int(dim.dataset[1], size = testSize)
    env <- cbind(knn$nn.index[idx.runs,1:(k0-1)], idx.runs)

    env.rand <- t(sapply(rep(dim.dataset[1],testSize),  sample.int, k0))

    #perform test
    p.val.test <- apply(env, 1, FUN = chi_batch_test, class.frequency, batch,  dof)
    p.val.test.null <- apply(env, 1, FUN = chi_batch_test, class.frequency, batch.shuff, dof)

    #summarise test results
    kBET.expected[i] <- sum(p.val.test.null < alpha) / length(p.val.test.null)
    kBET.observed[i] <- sum(p.val.test < alpha) / length(p.val.test)

    #compute significance
    kBET.signif[i] <- 1 - ptnorm(kBET.observed[i],
                                 mu=kBET.expected[i],
                                 sd=sqrt(kBET.expected[i]*(1-kBET.expected[i])/testSize),
                                 alpha=alpha)

    #assign results to result table
    rejection$results$tested[idx.runs] <- 1
    rejection$results$kBET.pvalue.test[idx.runs] <- p.val.test
    rejection$results$kBET.pvalue.null[idx.runs] <- p.val.test.null


    #compute likelihood-ratio test (approximation for multinomial exact test)
    p.val.test.lrt <- apply(env, 1, FUN = lrt_approximation, class.frequency, batch,  dof)
    p.val.test.lrt.null <- apply(env, 1, FUN = lrt_approximation, class.frequency, batch.shuff, dof)

    lrt.expected[i] <- sum(p.val.test.lrt.null < alpha) / length(p.val.test.lrt.null)
    lrt.observed[i] <- sum(p.val.test.lrt < alpha) / length(p.val.test.lrt)

    lrt.signif[i] <- 1 - ptnorm(lrt.observed[i],
                                mu=lrt.expected[i],
                                sd=sqrt(lrt.expected[i]*(1-lrt.expected[i])/testSize),
                                alpha=alpha)

    rejection$results$lrt.pvalue.test[idx.runs] <- p.val.test.lrt
    rejection$results$lrt.pvalue.null[idx.runs] <- p.val.test.lrt.null


    #exact result test consume a fairly high amount of computation time, as the exact test computes the probability of ALL
    #possible configurations (under the assumption that all batches are large enough to 'imitate' sampling with replacement)
    #For example: k0=33 and dof=5 yields 501942 possible choices and a computation time of several seconds (on a 16GB RAM machine)
    if (exists(x='exact.observed')){

      p.val.test.exact <- apply(env, 1, multiNom,class.frequency$freq, batch)

      p.val.test.exact.null <- apply(env, 1, multiNom, class.frequency$freq, batch.shuff)

      exact.expected[i] <- sum(p.val.test.exact.null<alpha)/testSize
      exact.observed[i] <- sum(p.val.test.exact<alpha)/testSize
      #compute the significance level for the number of rejected data points
      exact.signif[i] <- 1 - ptnorm(exact.observed[i],
                                    mu=exact.expected[i],
                                    sd=sqrt(exact.expected[i]*(1-exact.expected[i])/testSize),
                                    alpha=alpha)

      #p-value distribution
      rejection$results$exact.pvalue.test[idx.runs] <- p.val.test.exact
      rejection$results$exact.pvalue.null[idx.runs] <- p.val.test.exact.null
    }

  }



  if (stats>1){
    #summarize chi2-results
    CI95 <- c(0.025,0.5,0.975)
    rejection$summary$kBET.expected <-  c(mean(kBET.expected) ,quantile(kBET.expected, CI95))
    rownames(rejection$summary) <- c('mean', '2.5%', '50%', '97.5%')
    rejection$summary$kBET.observed <-  c(mean(kBET.observed) ,quantile(kBET.observed, CI95))
    rejection$summary$kBET.signif <- c(mean(kBET.signif) ,quantile(kBET.signif, CI95))
    #summarize lrt-results
    rejection$summary$lrt.expected <-  c(mean(lrt.expected) ,quantile(lrt.expected, CI95))

    rejection$summary$lrt.observed <-  c(mean(lrt.observed) ,quantile(lrt.observed, CI95))
    rejection$summary$lrt.signif <- c(mean(lrt.signif) ,quantile(lrt.signif, CI95))
    #names(rejection$lrt.observed)[1] <-'mean'
    #summarize exact test results
    if (exists(x='exact.observed')){
      rejection$summary$exact.expected <-  c(mean(exact.expected) ,quantile(exact.expected, CI95))
      rejection$summary$exact.observed <-  c(mean(exact.observed) ,quantile(exact.observed, CI95))
      rejection$summary$exact.signif <- c(mean(exact.signif) ,quantile(exact.signif, CI95))
    }

    if (stats<10){
      cat('Warning: The quantile computation for ')
      cat(paste0(stats))
      cat(' subset results is not meaningful.')
    }

    if(plot==TRUE & exists(x='exact.observed')){
      plot.data <- data.frame(class=rep(c('kBET', 'kBET (random)', 'lrt', 'lrt (random)', 'exact', 'exact (random)'), each=stats),
                              data =  c(kBET.observed, kBET.expected, lrt.observed, lrt.expected, exact.observed, exact.expected))
      g <-ggplot(plot.data, aes(class, data)) + geom_boxplot() +
        theme_bw() + labs(x='Test', y='Rejection rate')  + theme(axis.text.x = element_text(angle=45, hjust = 1))
      print(g)
    }
    if(plot==TRUE & !exists(x='exact.observed')){
      plot.data <- data.frame(class=rep(c('kBET', 'kBET (random)', 'lrt', 'lrt (random)'), each=stats),
                              data =  c(kBET.observed, kBET.expected, lrt.observed, lrt.expected))
      g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() +
        theme_bw() + labs(x='Test', y='Rejection rate') +
        scale_y_continuous(limits=c(0,1))  + theme(axis.text.x = element_text(angle=45, hjust = 1))
      print(g)
    }

  }else{ #i.e. no stats
    rejection$summary$kBET.expected <- kBET.expected
    rejection$summary$kBET.observed <- kBET.observed
    rejection$summary$kBET.signif <- kBET.signif

    if (addTest==TRUE){
      rejection$summary$lrt.expected <- lrt.expected
      rejection$summary$lrt.observed <- lrt.observed
      rejection$summary$lrt.signif <- lrt.signif
      if (exists(x='exact.observed')){
        rejection$summary$exact.expected <-  exact.expected
        rejection$summary$exact.observed <-  exact.observed
        rejection$summary$exact.signif <- exact.signif
      }
    }

  }
  }else{ #kBET only
    #initialize result list
    rejection <- list()
    rejection$summary <- data.frame(kBET.expected = numeric(4),
                                    kBET.observed = numeric(4),
                                    kBET.signif = numeric(4))
    kBET.expected <- numeric(stats)
    kBET.observed <- numeric(stats)
    kBET.signif <- numeric(stats)
    rejection$results   <- data.frame(tested = numeric(dim.dataset[1]),
                                      kBET.pvalue.test = rep(0,dim.dataset[1]),
                                      kBET.pvalue.null = rep(0, dim.dataset[1]))

      for (i in 1:stats){
      # choose a random sample from dataset (rows: samples, columns: parameters)
      idx.runs <- sample.int(dim.dataset[1], size = testSize)
      env <- cbind(knn$nn.index[idx.runs,1:(k0-1)], idx.runs)

      env.rand <- t(sapply(rep(dim.dataset[1],testSize),  sample.int, k0))

      #perform test
      p.val.test <- apply(env, 1, FUN = chi_batch_test, class.frequency, batch,  dof)
      p.val.test.null <- apply(env, 1, FUN = chi_batch_test, class.frequency, batch.shuff, dof)

      #summarise test results
      kBET.expected[i] <- sum(p.val.test.null < alpha) / length(p.val.test.null)
      kBET.observed[i] <- sum(p.val.test < alpha) / length(p.val.test)

      #compute significance
      kBET.signif[i] <- 1 - ptnorm(kBET.observed[i],
                                   mu=kBET.expected[i],
                                   sd=sqrt(kBET.expected[i]*(1-kBET.expected[i])/testSize),
                                   alpha=alpha)
      #assign results to result table
      rejection$results$tested[idx.runs] <- 1
      rejection$results$kBET.pvalue.test[idx.runs] <- p.val.test
      rejection$results$kBET.pvalue.null[idx.runs] <- p.val.test.null
    }

    if (stats>1){
      #summarize chi2-results
      CI95 <- c(0.025,0.5,0.975)
      rejection$summary$kBET.expected <-  c(mean(kBET.expected) ,quantile(kBET.expected, CI95))
      rownames(rejection$summary) <- c('mean', '2.5%', '50%', '97.5%')
      rejection$summary$kBET.observed <-  c(mean(kBET.observed) ,quantile(kBET.observed, CI95))
      rejection$summary$kBET.signif <- c(mean(kBET.signif) ,quantile(kBET.signif, CI95))

      #return also stats
      rejection$stats$kBET.expected <- kBET.expected
      rejection$stats$kBET.observed <- kBET.observed
      rejection$stats$kBET.signif <- kBET.signif

      if (stats<10){
        cat('Warning: The quantile computation for ')
        cat(paste0(stats))
        cat(' subset results is not meaningful.')
      }
      if(plot==TRUE){
        plot.data <- data.frame(class=rep(c('observed(kBET)', 'expected(random)'), each=stats), data =  c( kBET.observed,kBET.expected))
        g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + theme_bw() + labs(x='Test', y='Rejection rate') +
              scale_y_continuous(limits=c(0,1))
        print(g)
      }
    }else{
      rejection$summary$kBET.expected <- kBET.expected[1]
      rejection$summary$kBET.observed <- kBET.observed[1]
      rejection$summary$kBET.signif <- kBET.signif[1]
    }
  }

  return(rejection)
}

