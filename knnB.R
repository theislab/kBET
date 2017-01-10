# This function uses a chi square test to evaluate the probability of a
# batch effect. 

# df: dataset (rows: samples, columns: parameters)
# batch: batch id for each cell
# k0: number of nearest neighbours to test on (neighbourhood size)
# knn: a set of nearest neighbours for each cell (optional)
# testSize: number of data points to test, default: 10% sample size
# heuristic: compute an optimal neighbourhood size k 
# stats: to create a statistics on batch estimates, evaluate 'stats' subsets
# alpha: significance level
# addTest: perform an LRT-approximation to the multinomial test AND a multinomial exact test (if appropriate)
# plot: if stats > 10, then a boxplot of the resulting rejection rates is created
knnB <- function(df, batch, k0=NULL,knn=NULL, testSize=NULL,heuristic=FALSE, stats=100, alpha=0.05, addTest = FALSE, verbose=TRUE, plot = TRUE){
 require('ggplot2')
  
  if (plot==TRUE & heuristic==TRUE){
    source('addalpha.R')
    colorset <- addalpha(brewer.pal(8, 'Set2'), 0.5)
    colorset.nt <- brewer.pal(8, 'Set2')
  }
  
  if(addTest==TRUE){
    #load adapted EMT package
    path.name <- "EMT_adapted/" #how to adjust this? EMT_adapted was added to the repository
    file.sources = list.files(path= path.name, pattern= "*.R")
    sapply(paste0(path.name,file.sources),source,.GlobalEnv)
  }
  
  
  #preliminaries:
  dof <- length(unique(batch))-1 #degrees of freedom
  frequencies <- table(batch)/length(batch)
  class.frequency <- data.frame(class = names(frequencies), 
                                freq = as.numeric(frequencies))
  dataset <- df
  dim.dataset <- dim(dataset)
  stopifnot(class(stats) == 'numeric', stats>0)
  
  if (is.null(k0) || k0>=dim.dataset[1]){
    #default environment size: quarter the size of the largest batch
    k0=floor(mean(class.frequency$freq)*dim.dataset[1]/4)
    cat('size of neighbourhood is set to ')
    cat(paste0(k0, '\n'))
  }
  # find KNNs
  if (is.null(knn)){
    require('FNN') #nearest neighbour search
    if (verbose) {
      cat('finding knns...')
      tic <- proc.time()
    }
    knn <- get.knn(dataset, k=k0, algorithm = 'cover_tree')
    if (verbose) {
      cat('...done. Time:\n')
      print(proc.time() - tic)
    }
  }
  #set number of tests 
  if (is.null(testSize) || dim.dataset[1]< testSize){
    test.frac <- 0.1
    testSize <- ceiling(dim.dataset[1]*test.frac)
    cat('number of neighbourhood tests is set to ')
    cat(paste0(testSize, '\n'))
  }
  
  if(heuristic==TRUE){
    require(TTR)
    test.k <- 30
    k <- round(seq(10, k0, length.out = test.k),0)
    myfun <- function(x,df,batch, knn){
      res <- knnB(df=df, batch=batch, k0=x, knn=knn, testSize=NULL, heuristic=FALSE, stats=10, alpha=0.05, addTest = FALSE, plot=FALSE)
      result <- res$summary
    }
    #btw, when we use sapply here in this context, it creates a list.
    heuristics <- sapply(k, myfun, df, batch, knn)
    data.heu <- matrix(unlist(heuristics['knnB.observed',]), ncol=test.k)
    #analysis of heuristics - take the average max of some moving average
    window.size  <- 5
    mean.reject <- apply(X=data.heu, MARGIN=1, FUN = SMA, window.size)
    pre.opt.k <- window.size-1 + apply(mean.reject[window.size: (length(k)),], 2, function(x) {max(which(x==max(x)))})
    opt.k <- floor(mean(pre.opt.k))
    #result
    k0 <- k[opt.k]
    if(verbose==TRUE){
      cat('The optimal neighbourhood size is determined.\n')
      cat('Size of neighbourhood is set to ')
      cat(paste0(k0, '\n'))
    }
    if (plot==TRUE){
      k.frac <- 100*k/dim(df)[1]
      data.null <- matrix(unlist(heuristics['knnB.expected',]), ncol=test.k)
      data.kdep <- data.frame(neighbourhood=rep(c(k.frac, rev(k.frac)), 2), 
                              result=c(data.heu[2,], rev(data.heu[4,]), #chi2 observed
                                       data.null[2,], rev(data.null[4,])), #chi2 expected
                              test= as.factor(rep(c('knnB','chi2 expected'), each=2*length(k)))) 
      data.kdep.points <- data.frame(neighbourhood=rep(k.frac,2), result= c(data.heu[1,], data.null[1,]),
                                     test= as.factor(rep(c('knnB','chi2 expected'), each=length(k))))
      
      p.heu <- ggplot(data.kdep, aes(x=neighbourhood, y=result)) +
        geom_polygon(data=data.kdep, mapping=aes(x=neighbourhood, y=result, group=test, fill=test)) +
        geom_line(data=data.kdep.points, mapping=aes(x=neighbourhood, y=result, color=test)) +
        scale_color_manual(values = c('knnB' = colorset[1], 'chi2 expected'= colorset[2]),
                           name = 'Test', breaks=c('knnB','chi2 expected'),
                           labels = c('knnB', 'random\nassignment')) +
        scale_fill_manual(values = c('knnB' = colorset[1], 'chi2 expected'= colorset[2]),
                          name = 'Test', breaks=c('knnB','chi2 expected'),
                          labels = c('knnB', 'random\nassignment')) +
        labs(x= 'Neighbourhood size (in % sample size)', y = 'Rejection rate')+ 
        geom_vline(aes(xintercept = k.frac[opt.k]),
                   linetype="dotted" , size=0.5)+
        theme_bw()
    }
    print(p.heu)
  }
  
  if(addTest==TRUE){
    #initialize result list  
  rejection <- list()
  rejection$summary <- data.frame(knnB.expected = numeric(4),
                                  knnB.observed = numeric(4),
                                  lrt.expected = numeric(4),
                                  lrt.observed = numeric(4)
  ) 
  knnB.expected <- numeric(stats)
  knnB.observed <- numeric(stats)
  knnB.signif <- numeric(stats)
  rejection$results   <- data.frame(tested = numeric(dim.dataset[1]),
                                    knnB.pvalue.test = rep(0,dim.dataset[1]),
                                    knnB.pvalue.null = rep(0, dim.dataset[1]),
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
    p.val.test.null <- apply(env.rand, 1, FUN = chi_batch_test, class.frequency, batch, dof) 
    
    #summarise test results
    knnB.expected[i] <- sum(p.val.test.null < alpha) / length(p.val.test.null)
    knnB.observed[i] <- sum(p.val.test < alpha) / length(p.val.test)
    
    #compute significance
    knnB.signif[i] <- 1 - ptnorm(knnB.observed[i], 
                                 mu=knnB.expected[i], 
                                 sd=sqrt(knnB.expected[i]*(1-knnB.expected[i])/testSize), 
                                 alpha=alpha)
    #assign results to result table
    rejection$results$tested[idx.runs] <- 1
    rejection$results$knnB.pvalue.test[idx.runs] <- p.val.test
    rejection$results$knnB.pvalue.null[idx.runs] <- p.val.test.null
    
      
    #compute likelihood-ratio test (approximation for multinomial exact test)
    p.val.test.lrt <- apply(env, 1, FUN = lrt_approximation, class.frequency, batch,  dof)
    p.val.test.lrt.null <- apply(env.rand, 1, FUN = lrt_approximation, class.frequency, batch, dof)
    
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
      
      p.val.test.exact <- apply(env, 1, multiNom, 
        class.frequency$freq, batch)
      
      p.val.test.exact.null <- apply(env.rand, 1, multiNom, class.frequency$freq, batch)
      
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
      if(plot==TRUE){
        plot.data <- data.frame(class=rep(c('knnB', 'knnB (random)', 'lrt', 'lrt (random)', 'exact', 'exact (random)'), each=stats), 
                                data =  c(knnB.observed, knnB.expected, lrt.observed, lrt.expected, exact.observed, exact.expected))
        ggplot(plot.data, aes(class, data)) + geom_boxplot() + theme_bw() + labs(x='Test', y='Rejection rate')
      }
      if(plot==TRUE & !exists(x='exact.observed')){
        plot.data <- data.frame(class=rep(c('knnB', 'knnB (random)', 'lrt', 'lrt (random)'), each=stats), 
                                data =  c(knnB.observed, knnB.expected, lrt.observed, lrt.expected))
        g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + theme_bw() + labs(x='Test', y='Rejection rate') + scale_y_continuous(limits=c(0,1))
        print(g)
      }
    }
    }
   
   
  
  if (stats>1){
    #summarize chi2-results
    CI95 <- c(0.025,0.5,0.975)
    rejection$summary$knnB.expected <-  c(mean(knnB.expected) ,quantile(knnB.expected, CI95))
    rownames(rejection$summary) <- c('mean', '2.5%', '50%', '97.5%')
    rejection$summary$knnB.observed <-  c(mean(knnB.observed) ,quantile(knnB.observed, CI95))
    rejection$summary$knnB.signif <- c(mean(knnB.signif) ,quantile(knnB.signif, CI95))
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
  }else{ #i.e. no stats
    rejection$summary$knnB.expected <- knnB.expected
    rejection$summary$knnB.observed <- knnB.observed
    rejection$summary$knnB.signif <- knnB.signif

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
  }else{ #knnB only
    #initialize result list
    rejection <- list()
    rejection$summary <- data.frame(knnB.expected = numeric(4),
                                    knnB.observed = numeric(4),
                                    knnB.signif = numeric(4)) 
    knnB.expected <- numeric(stats)
    knnB.observed <- numeric(stats)
    knnB.signif <- numeric(stats)
    rejection$results   <- data.frame(tested = numeric(dim.dataset[1]),
                                      knnB.pvalue.test = rep(0,dim.dataset[1]),
                                      knnB.pvalue.null = rep(0, dim.dataset[1]))
  
      for (i in 1:stats){
      # choose a random sample from dataset (rows: samples, columns: parameters)
      idx.runs <- sample.int(dim.dataset[1], size = testSize)
      env <- cbind(knn$nn.index[idx.runs,1:(k0-1)], idx.runs)
      
      env.rand <- t(sapply(rep(dim.dataset[1],testSize),  sample.int, k0))
      
      #perform test
      p.val.test <- apply(env, 1, FUN = chi_batch_test, class.frequency, batch,  dof)
      p.val.test.null <- apply(env.rand, 1, FUN = chi_batch_test, class.frequency, batch, dof) 
      
      #summarise test results
      knnB.expected[i] <- sum(p.val.test.null < alpha) / length(p.val.test.null)
      knnB.observed[i] <- sum(p.val.test < alpha) / length(p.val.test)
      
      #compute significance
      knnB.signif[i] <- 1 - ptnorm(knnB.observed[i], 
                                   mu=knnB.expected[i], 
                                   sd=sqrt(knnB.expected[i]*(1-knnB.expected[i])/testSize), 
                                   alpha=alpha)
      #assign results to result table
      rejection$results$tested[idx.runs] <- 1
      rejection$results$knnB.pvalue.test[idx.runs] <- p.val.test
      rejection$results$knnB.pvalue.null[idx.runs] <- p.val.test.null
    }
 
    if (stats>1){
      #summarize chi2-results
      CI95 <- c(0.025,0.5,0.975)
      rejection$summary$knnB.expected <-  c(mean(knnB.expected) ,quantile(knnB.expected, CI95))
      rownames(rejection$summary) <- c('mean', '2.5%', '50%', '97.5%')
      rejection$summary$knnB.observed <-  c(mean(knnB.observed) ,quantile(knnB.observed, CI95))
      rejection$summary$knnB.signif <- c(mean(knnB.signif) ,quantile(knnB.signif, CI95))
      
      if (stats<10){
        cat('Warning: The quantile computation for ')
        cat(paste0(stats))
        cat(' subset results is not meaningful.')
      }
      if(plot==TRUE){
        plot.data <- data.frame(class=rep(c('knnB', 'random assigment'), each=stats), data =  c( knnB.observed,knnB.expected))
        g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + theme_bw() + labs(x='Test', y='Rejection rate') + 
              scale_y_continuous(limits=c(0,1))
        print(g)
      }
    }else{
      rejection$summary$knnB.expected <- knnB.expected[1]
      rejection$summary$knnB.observed <- knnB.observed[1]
      rejection$summary$knnB.signif <- knnB.signif[1]
    }
  }

  return(rejection)
}

chi_batch_test <- function(knn.set, class.freq, batch, df)
{
  #knn.set: indices of nearest neighbours
  #empirical frequencies in nn-environment (sample 1)
  freq.env <- table(batch[knn.set])
  full.classes <- rep(0, length(class.freq$class))
  full.classes[ class.freq$class %in% names(freq.env)] <- freq.env
  exp.freqs <- class.freq$freq*length(knn.set)
  #compute chi-square test statistics
  chi.sq.value <- sum((full.classes - exp.freqs)^2/exp.freqs)
  result<- 1- pchisq(chi.sq.value, df) #p-value for the result
  return(result)
}

lrt_approximation <- function(knn.set, class.freq, batch, df)
{
  #knn.set: indices of nearest neighbours
  #empirical frequencies in nn-environment (sample 1)
  obs.env <- table(batch[knn.set]) #observed realisations of each category
  freq.env <- obs.env/sum(obs.env) #observed 'probabilities'
  full.classes <- rep(0, length(class.freq$class))
  obs.classes <- class.freq$class %in% names(freq.env)
  #for stability issues (to avoid the secret division by 0): introduce another alternative model where the observed probability
  #is either the empirical frequency or 1/(sample size) at minimum
  if (length(full.classes) > sum(obs.classes)){
    dummy.count <- length(full.classes) -sum(obs.classes)
    full.classes[obs.classes] <- obs.env/(sum(obs.env)+ dummy.count)
    pmin <- 1/(sum(obs.env)+ dummy.count)
    full.classes[!obs.classes] <- pmin
  }else{
    full.classes[ obs.classes] <- freq.env
  }
  exp.freqs <- class.freq$freq #expected 'probabilities'
  #compute likelihood ratio of null and alternative hypothesis, test statistics converges to chi-square distribution
  full.obs <- rep(0, length(class.freq$class))
  full.obs[obs.classes] <- obs.env
  
  lrt.value <- -2*sum(full.obs * log(exp.freqs/full.classes))
  
  result<- 1- pchisq(lrt.value, df) #p-value for the result
  return(result)
}
#truncated normal distribution distribution function
ptnorm <- function(x,mu,sd, a=0, b=1, alpha=0.05){
  #this is the cumulative density of the truncated normal distribution
  #x ~ N(mu, sd^2), but we condition on a <= x <= b
  if(a>b){
    warning("Lower and upper bound are interchanged.")
    tmp <- a
    a <- b
    b <- tmp 
  }
  
  if(sd<=0) {
    warning("Standard deviation must be positive.")
    if (alpha<=0)
    {
      stop("False positive rate alpha must be positive.")
    }
    sd <- alpha
  } 
  if (x<a | x>b){
    warning("x out of bounds.")
    cdf <- as.numeric(x>a)
  }else{
    alp <- pnorm((a-mu)/sd)
    bet <- pnorm((b-mu)/sd)
    zet <- pnorm((x-mu)/sd)
    cdf <- (zet-alp)/(bet-alp)
  }
  return(cdf)
}
#wrapper for the multinomial exact test function
multiNom <- function(x, y, z) {
  z.f <- factor(z)
  tmp <- multinomial.test(as.numeric(table(z.f[x])),y)
  return(tmp$p.value)}
