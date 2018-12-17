#' kBET - k-nearest neighbour batch effect test
#'
#' @description \code{kBET} runs a chi square test to evaluate
#' the probability of a batch effect.
#'
#' @param df dataset (rows: samples, columns: features)
#' @param batch batch id for each cell or a data frame with
#' both condition and replicates
#' @param k0 number of nearest neighbours to test on (neighbourhood size)
#' @param knn a set of nearest neighbours for each cell (optional)
#' @param testSize number of data points to test,
#' (10 percent sample size default, but at least 25)
#' @param do.pca perform a pca prior to knn search? (defaults to TRUE)
#' @param dim.pca if do.pca=TRUE, choose the number of dimensions
#' to consider (defaults to 50)
#' @param heuristic compute an optimal neighbourhood size k
#' (defaults to TRUE)
#' @param n_repeat to create a statistics on batch estimates,
#' evaluate 'n_repeat' subsets
#' @param alpha significance level
#' @param adapt In some cases, a number of cells do not contribute
#' to any neighbourhood
#' and this may cause an imbalance in observed and expected batch
#' label frequencies.
#' Frequencies will be adapted if adapt=TRUE (default).
#' @param addTest perform an LRT-approximation to the multinomial
#' test AND a multinomial exact test (if appropriate)
#' @param plot if stats > 10, then a boxplot of the resulting
#' rejection rates is created
#' @param verbose displays stages of current computation (defaults to FALSE)
#' @return list object
#'    \enumerate{
#'    \item \code{summary} - a rejection rate for the data,
#'         an expected rejection rate for random
#'         labeling and the significance for the observed result
#'    \item \code{results} - detailed list for each tested cells;
#'         p-values for expected and observed label distribution
#'    \item \code{average.pval} - significance level over the averaged
#'         batch label distribution in all neighbourhoods
#'    \item \code{stats} - extended test summary for every sample
#'    \item \code{params} - list of input parameters and adapted parameters,
#'    respectively
#'    \item \code{outsider} - only shown if \code{adapt=TRUE}. List of samples
#'         without mutual nearest neighbour: \itemize{
#'     \item \code{index} - index of each outsider sample)
#'     \item \code{categories} - tabularised labels of outsiders
#'     \item \code{p.val} - Significance level of outsider batch label distribution
#'         vs expected frequencies.
#'     If the significance level is lower than \code{alpha},
#'     expected frequencies will be adapted}
#'    }
#' @examples
#'     batch <- rep(seq_len(10),each=20)
#'     data <- matrix(rpois(n = 50000, lambda = 10)*rbinom(50000,1,prob=0.5), ncol=200)
#'
#'     batch.estimate <- kBET(data,batch)

#' @importFrom FNN get.knn
#' @import ggplot2
#' @importFrom stats quantile pchisq
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils data
#' @include bisect.R
#' @include kBET-utils.R
#' @name kBET
#' @export
kBET <- function(df, batch, k0=NULL, knn = NULL,
                 testSize=NULL, do.pca = TRUE, dim.pca = 50,
                 heuristic=TRUE, n_repeat = 100,
                 alpha=0.05, addTest = FALSE,
                 verbose=FALSE, plot = TRUE, adapt=TRUE) {

  #create a subsetting mode:
  #if (is.data.frame(batch) && dim(batch)[2]>1) {
  #  cat('Complex study design detected.\n')
  #  design.names <- colnames(batch)
  #  cat(paste0('Subset by ', design.names[[2]], '\n'))
  #  bio <- unlist(unique(batch[[design.names[2]]]))
  #drop subset batch
  #  new.batch <- base::subset(batch,
  # batch[[2]]== bio[1])[,!design.names %in% design.names[[2]]]
  #  sapply(df[batch[[2]]==bio[1],], kBET, new.batch,
  # k0, knn, testSize, heuristic,n_repeat, alpha, addTest, verbose)

  #}


  #preliminaries:
  dof <- length(unique(batch)) - 1 #degrees of freedom
  if (is.factor(batch)) {
    batch <- droplevels(batch)
  }

  frequencies <- table(batch)/length(batch)
  #get 3 different permutations of the batch label
  batch.shuff <- replicate(3, batch[sample.int(length(batch))])


  class.frequency <- data.frame(class = names(frequencies),
                                freq = as.numeric(frequencies))
  dataset <- df
  dim.dataset <- dim(dataset)
  #check the feasibility of data input
  if (dim.dataset[1] != length(batch) && dim.dataset[2] != length(batch)) {
    stop("Input matrix and batch information do not match. Execution halted.")
  }

  if (dim.dataset[2] == length(batch) && dim.dataset[1] != length(batch)) {
    if (verbose) {
      cat('Input matrix has samples as columns. kBET needs samples as rows.
        Transposing...\n')
    }
    dataset <- t(dataset)
    dim.dataset <- dim(dataset)
  }


  stopifnot(class(n_repeat) == 'numeric', n_repeat > 0)

  do_heuristic <- FALSE
  if (is.null(k0) || k0 >= dim.dataset[1]) {
    do_heuristic <- heuristic
    if (!heuristic) {
      #default environment size: quarter the size of the largest batch
      k0 <- floor(mean(class.frequency$freq)*dim.dataset[1]/4)
    } else {
      #default environment size: three quarter the size of the largest batch
      k0 <- floor(mean(class.frequency$freq)*dim.dataset[1]*0.75)
      if (k0 < 10) {
        stop("Your dataset has too few samples to run a heuristic.\n
             Please assign k0 and set heuristic=FALSE.")
      }
    }
    if (verbose) {
      cat('Initial neighbourhood size is set to ')
      cat(paste0(k0, '.\n'))
    }
  }
  # find KNNs
  if (is.null(knn)) {
    if (!do.pca) {
      if (verbose) {
        cat('finding knns...')
        tic <- proc.time()
      }
      knn <- get.knn(dataset, k = k0, algorithm = 'cover_tree')
    } else {
      dim.comp <- min(dim.pca, dim.dataset[2])
      if (verbose) {
        cat('reducing dimensions with svd first...\n')
      }
      data.pca <- svd(x = dataset, nu = dim.comp, nv = 0)
      if (verbose) {
        cat('finding knns...')
        tic <- proc.time()
      }
      knn <- get.knn(data.pca$u,  k = k0, algorithm = 'cover_tree')
    }
    if (verbose) {
      cat('done. Time:\n')
      print(proc.time() - tic)
    }
  }

  #set number of tests
  if (is.null(testSize) || (floor(testSize) < 1 || dim.dataset[1] < testSize)) {
    test.frac <- 0.1
    testSize <- ceiling(dim.dataset[1]*test.frac)
    if (testSize < 25 && dim.dataset[1] > 25) {
      testSize <- 25
    }
    if (verbose) {
      cat('Number of kBET tests is set to ')
      cat(paste0(testSize, '.\n'))
    }
  }
  #decide to adapt general frequencies
  if (adapt) {
    # idx.run <- sample.int(dim.dataset[1], size = min(2*testSize, dim.dataset[1]))
    outsider <- which(!(seq_len(dim.dataset[1]) %in% knn$nn.index[, seq_len(k0 - 1)]))
    is.imbalanced <- FALSE #initialisation
    p.out <- 1
    #avoid unwanted things happen if length(outsider) == 0
    if (length(outsider) > 0) {
      p.out <- chi_batch_test(outsider, class.frequency, batch,  dof)
      is.imbalanced <- p.out < alpha
      if (is.imbalanced) {
        new.frequencies <- table(batch[-outsider])/length(batch[-outsider])
        new.class.frequency <- data.frame(class = names(new.frequencies),
                                          freq = as.numeric(new.frequencies))
        if (verbose) {
          cat(paste0('There are ', length(outsider), ' cells (',
                     round(length(outsider)/length(batch)*100,3),
                     '%) that do not appear in any neighbourhood.\n',
                     'The expected frequencies for each category have
                     been adapted.\n',
                     'Cell indexes are saved to result list.\n'))
        }
      } else {
        if (verbose) {
          cat(paste0('No outsiders found.'))
        }
      }
    }

  }

  if (do_heuristic) {
    #btw, when we bisect here that returns some interval with
    #the optimal neihbourhood size
    if (verbose) {
      cat('Determining optimal neighbourhood size ...')
    }
    opt.k <- bisect(scan_nb, bounds = c(10, k0), known = NULL, dataset, batch, knn)
    #result
    if (length(opt.k) > 1) {
      k0 <- opt.k[2]
      if (verbose) {
        cat(paste0('done.\nNew size of neighbourhood is set to ', k0, '.\n'))
      }
    } else {
      if (verbose) {
        cat(paste0('done.\nHeuristic did not change the
                   neighbourhood.\n If results appear inconclusive,
                   change k0=', k0, '.\n'))
      }
    }
  }

  #initialise result list
  rejection <- list()
  rejection$summary <- data.frame(kBET.expected = numeric(4),
                                  kBET.observed = numeric(4),
                                  kBET.signif = numeric(4))

  rejection$results   <- data.frame(tested = numeric(dim.dataset[1]),
                                    kBET.pvalue.test = rep(0,dim.dataset[1]),
                                    kBET.pvalue.null = rep(0, dim.dataset[1]))
  #get average residual score
  env <- as.vector(cbind(knn$nn.index[, seq_len(k0 - 1)], seq_len(dim.dataset[1])))
  cf <- if (adapt && is.imbalanced) new.class.frequency else class.frequency
  rejection$average.pval <- 1 - pchisq(k0 * residual_score_batch(env, cf, batch), dof)



  #initialise intermediates
  kBET.expected <- numeric(n_repeat)
  kBET.observed <- numeric(n_repeat)
  kBET.signif <- numeric(n_repeat)


  if (addTest) {
    #initialize result list
    rejection$summary$lrt.expected <- numeric(4)
    rejection$summary$lrt.observed <- numeric(4)

    rejection$results$lrt.pvalue.test <- rep(0,dim.dataset[1])
    rejection$results$lrt.pvalue.null <- rep(0, dim.dataset[1])

    lrt.expected <- numeric(n_repeat)
    lrt.observed <- numeric(n_repeat)
    lrt.signif <- numeric(n_repeat)
    #decide to perform exact test or not
    if (choose(k0 + dof,dof) < 5e5 && k0 <= min(table(batch))) {
      exact.expected <- numeric(n_repeat)
      exact.observed <- numeric(n_repeat)
      exact.signif <- numeric(n_repeat)

      rejection$summary$exact.expected <- numeric(4)
      rejection$summary$exact.observed <- numeric(4)
      rejection$results$exact.pvalue.test <- rep(0, dim.dataset[1])
      rejection$results$exact.pvalue.null <- rep(0, dim.dataset[1])
    }

    for (i in seq_len(n_repeat)) {
      # choose a random sample from dataset (rows: samples, columns: features)
      idx.runs <- sample.int(dim.dataset[1], size = testSize)
      env <- cbind(knn$nn.index[idx.runs, seq_len(k0 - 1)], idx.runs)
      #env.rand <- t(sapply(rep(dim.dataset[1],testSize),  sample.int, k0))

      #perform test
      if (adapt && is.imbalanced) {
        p.val.test <- apply(env, 1, FUN = chi_batch_test,
                            new.class.frequency, batch,  dof)
      } else {
        p.val.test <- apply(env, 1, FUN = chi_batch_test,
                            class.frequency, batch,  dof)
      }

      is.rejected <- p.val.test < alpha



      #p.val.test <- apply(env, 1, FUN = chi_batch_test, class.frequency,
      #batch,  dof)
      p.val.test.null <-  apply(apply(batch.shuff, 2,
                                      function(x, freq, dof, envir) {
                                        apply(envir, 1, FUN = chi_batch_test, freq, x, dof)},
                                      class.frequency, dof, env), 1, mean)
      #p.val.test.null <- apply(env.rand, 1, FUN = chi_batch_test,
      #class.frequency, batch, dof)

      #summarise test results
      kBET.expected[i] <- sum(p.val.test.null < alpha) / length(p.val.test.null)
      kBET.observed[i] <- sum(is.rejected) / length(p.val.test)

      #compute significance
      kBET.signif[i] <-
        1 - ptnorm(kBET.observed[i],
                   mu = kBET.expected[i],
                   sd = sqrt(kBET.expected[i] * (1 - kBET.expected[i]) / testSize),
                   alpha = alpha)

      #assign results to result table
      rejection$results$tested[idx.runs] <- 1
      rejection$results$kBET.pvalue.test[idx.runs] <- p.val.test
      rejection$results$kBET.pvalue.null[idx.runs] <- p.val.test.null


      #compute likelihood-ratio test (approximation for multinomial exact test)
      cf <- if (adapt && is.imbalanced) new.class.frequency else class.frequency
      p.val.test.lrt <- apply(env, 1, FUN = lrt_approximation, new.class.frequency, batch, dof)
      p.val.test.lrt.null <- apply(apply(batch.shuff, 2,
                                         function(x, freq, dof, envir) {
                                           apply(envir, 1, FUN = lrt_approximation, freq, x, dof)},
                                         class.frequency, dof, env), 1, mean)

      lrt.expected[i] <-
        sum(p.val.test.lrt.null < alpha) / length(p.val.test.lrt.null)
      lrt.observed[i] <-
        sum(p.val.test.lrt < alpha) / length(p.val.test.lrt)

      lrt.signif[i] <-
        1 - ptnorm(lrt.observed[i],
                   mu = lrt.expected[i],
                   sd = sqrt(lrt.expected[i] * (1 - lrt.expected[i]) / testSize),
                   alpha = alpha)

      rejection$results$lrt.pvalue.test[idx.runs] <- p.val.test.lrt
      rejection$results$lrt.pvalue.null[idx.runs] <- p.val.test.lrt.null


      #exact result test consume a fairly high amount of computation time,
      #as the exact test computes the probability of ALL
      #possible configurations (under the assumption that all batches
      #are large enough to 'imitate' sampling with replacement)
      #For example: k0=33 and dof=5 yields 501942 possible
      #choices and a computation time of several seconds (on a 16GB RAM machine)
      if (exists(x = 'exact.observed')) {
        if (adapt && is.imbalanced) {
          p.val.test.exact <- apply(env, 1, multiNom,
                                    new.class.frequency$freq, batch)
        } else {
          p.val.test.exact <- apply(env, 1, multiNom,
                                    class.frequency$freq, batch)
        }
        p.val.test.exact.null <-  apply(apply(batch.shuff, 2,
                                              function(x, freq, envir) {
                                                apply(envir, 1, FUN = multiNom, freq, x)},
                                              class.frequency$freq, env), 1, mean)
        # apply(env, 1, multiNom, class.frequency$freq, batch.shuff)

        exact.expected[i] <- sum(p.val.test.exact.null < alpha)/testSize
        exact.observed[i] <- sum(p.val.test.exact < alpha)/testSize
        #compute the significance level for the number of rejected data points
        exact.signif[i] <-
          1 - ptnorm(exact.observed[i],
                     mu = exact.expected[i],
                     sd = sqrt(exact.expected[i] * (1 - exact.expected[i]) / testSize),
                     alpha = alpha)
        #p-value distribution
        rejection$results$exact.pvalue.test[idx.runs] <- p.val.test.exact
        rejection$results$exact.pvalue.null[idx.runs] <- p.val.test.exact.null
      }
    }



    if (n_repeat > 1) {
      #summarize chi2-results
      CI95 <- c(0.025,0.5,0.975)
      rejection$summary$kBET.expected <-  c(mean(kBET.expected) ,
                                            quantile(kBET.expected, CI95))
      rownames(rejection$summary) <- c('mean', '2.5%', '50%', '97.5%')
      rejection$summary$kBET.observed <-  c(mean(kBET.observed) ,
                                            quantile(kBET.observed, CI95))
      rejection$summary$kBET.signif <- c(mean(kBET.signif) ,
                                         quantile(kBET.signif, CI95))
      #summarize lrt-results
      rejection$summary$lrt.expected <-  c(mean(lrt.expected) ,
                                           quantile(lrt.expected, CI95))
      rejection$summary$lrt.observed <-  c(mean(lrt.observed) ,
                                           quantile(lrt.observed, CI95))
      rejection$summary$lrt.signif <- c(mean(lrt.signif) ,
                                        quantile(lrt.signif, CI95))
      #summarize exact test results
      if (exists(x = 'exact.observed')) {
        rejection$summary$exact.expected <-  c(mean(exact.expected) ,
                                               quantile(exact.expected, CI95))
        rejection$summary$exact.observed <-  c(mean(exact.observed) ,
                                               quantile(exact.observed, CI95))
        rejection$summary$exact.signif <- c(mean(exact.signif) ,
                                            quantile(exact.signif, CI95))
      }

      if (n_repeat < 10) {
        cat('Warning: The quantile computation for ')
        cat(paste0(n_repeat))
        cat(' subset results is not meaningful.')
      }

      if (plot && exists(x = 'exact.observed')) {
        plot.data <- data.frame(class = rep(c('kBET', 'kBET (random)',
                                            'lrt', 'lrt (random)',
                                            'exact', 'exact (random)'),
                                          each = n_repeat),
                                data =  c(kBET.observed, kBET.expected,
                                          lrt.observed, lrt.expected,
                                          exact.observed, exact.expected))
        g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() +
          theme_bw() + labs(x = 'Test', y = 'Rejection rate')  +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(g)
      }
      if (plot && !exists(x = 'exact.observed')) {
        plot.data <- data.frame(class = rep(c('kBET', 'kBET (random)',
                                            'lrt', 'lrt (random)'),
                                          each = n_repeat),
                                data =  c(kBET.observed, kBET.expected,
                                          lrt.observed, lrt.expected))
        g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() +
          theme_bw() + labs(x = 'Test', y  = 'Rejection rate') +
          scale_y_continuous(limits = c(0,1))  +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(g)
      }

    } else {#i.e. no n_repeat
      rejection$summary$kBET.expected <- kBET.expected
      rejection$summary$kBET.observed <- kBET.observed
      rejection$summary$kBET.signif <- kBET.signif

      if (addTest) {
        rejection$summary$lrt.expected <- lrt.expected
        rejection$summary$lrt.observed <- lrt.observed
        rejection$summary$lrt.signif <- lrt.signif
        if (exists(x = 'exact.observed')) {
          rejection$summary$exact.expected <-  exact.expected
          rejection$summary$exact.observed <-  exact.observed
          rejection$summary$exact.signif <- exact.signif
        }
      }

    }
  } else {#kBET only
    for (i in seq_len(n_repeat)) {
      # choose a random sample from dataset
      #(rows: samples, columns: parameters)
      idx.runs <- sample.int(dim.dataset[1], size = testSize)
      env <- cbind(knn$nn.index[idx.runs,seq_len(k0 - 1)], idx.runs)

      #perform test
      if (adapt && is.imbalanced) {
        p.val.test <- apply(env, 1, FUN = chi_batch_test,
                            new.class.frequency, batch,  dof)
      } else {
        p.val.test <- apply(env, 1, FUN = chi_batch_test,
                            class.frequency, batch,  dof)
      }

      #print(dim(env))
      is.rejected <- p.val.test < alpha

      p.val.test.null <- apply(batch.shuff, 2,
                               function(x, freq, dof, envir) {
                                 apply(envir, 1, FUN = chi_batch_test, freq, x, dof)},
                               class.frequency, dof, env)
      # p.val.test.null <- apply(env, 1, FUN = chi_batch_test,
      #class.frequency, batch.shuff, dof)

      #summarise test results
      #kBET.expected[i] <- sum(p.val.test.null < alpha) / length(p.val.test.null)
      kBET.expected[i] <- mean(apply(p.val.test.null, 2,
                                     function(x,alpha) {
                                       sum(x < alpha) / length(x)},
                                     alpha))

      kBET.observed[i] <- sum(is.rejected) / length(p.val.test)

      #compute significance
      kBET.signif[i] <-
        1 - ptnorm(kBET.observed[i],
                   mu = kBET.expected[i],
                   sd = sqrt(kBET.expected[i] * (1 - kBET.expected[i]) / testSize),
                   alpha = alpha)
      #assign results to result table
      rejection$results$tested[idx.runs] <- 1
      rejection$results$kBET.pvalue.test[idx.runs] <- p.val.test
      rejection$results$kBET.pvalue.null[idx.runs] <- rowMeans(p.val.test.null)
    }

    if (n_repeat > 1) {
      #summarize chi2-results
      CI95 <- c(0.025,0.5,0.975)
      rejection$summary$kBET.expected <- c(mean(kBET.expected), quantile(kBET.expected, CI95))
      rownames(rejection$summary) <- c('mean', '2.5%', '50%', '97.5%')
      rejection$summary$kBET.observed <- c(mean(kBET.observed), quantile(kBET.observed, CI95))
      rejection$summary$kBET.signif <- c(mean(kBET.signif), quantile(kBET.signif, CI95))

      #return also n_repeat
      rejection$stats$kBET.expected <- kBET.expected
      rejection$stats$kBET.observed <- kBET.observed
      rejection$stats$kBET.signif <- kBET.signif

      if (n_repeat < 10) {
        cat('Warning: The quantile computation for ')
        cat(paste0(n_repeat))
        cat(' subset results is not meaningful.')
      }
      if (plot) {
        plot.data <- data.frame(class = rep(c('observed(kBET)',
                                            'expected(random)'),
                                          each = n_repeat),
                                data =  c( kBET.observed,kBET.expected))
        g <- ggplot(plot.data, aes(class, data)) +
          geom_boxplot() + theme_bw() +
          labs(x = 'Test', y = 'Rejection rate') +
          scale_y_continuous(limits = c(0,1))
        print(g)
      }
    } else {
      rejection$summary$kBET.expected <- kBET.expected[1]
      rejection$summary$kBET.observed <- kBET.observed[1]
      rejection$summary$kBET.signif <- kBET.signif[1]
    }
  }
  #collect parameters
  rejection$params <- list()
  rejection$params$k0 <- k0
  rejection$params$testSize <- testSize
  rejection$params$do.pca <- do.pca
  rejection$params$dim.pca <- dim.pca
  rejection$params$heuristic <- heuristic
  rejection$params$n_repeat <- n_repeat
  rejection$params$alpha <- alpha
  rejection$params$addTest <- addTest
  rejection$params$verbose <- verbose
  rejection$params$plot <- plot

  #add outsiders
  if (adapt) {
    rejection$outsider <- list()
    rejection$outsider$index <- outsider
    rejection$outsider$categories <- table(batch[outsider])
    rejection$outsider$p.val <- p.out
  }
  rejection
}

