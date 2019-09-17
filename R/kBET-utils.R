#' @importFrom stats pchisq pnorm
#a wrapper for kBET to fix a neighbourhood size
scan_nb <- function(x,df,batch, knn) {
    res <- kBET(df = df, batch = batch, k0 = x, knn = knn, testSize = NULL,
                heuristic = FALSE, n_repeat = 10, alpha = 0.05,
                addTest = FALSE, plot = FALSE, verbose = FALSE, adapt = FALSE)
    result <- res$summary
    result$kBET.observed[1]
}

#the residual score function of kBET
residual_score_batch <- function(knn.set, class.freq, batch) {
  #knn.set: indices of nearest neighbours
  #empirical frequencies in nn-environment (sample 1)
  #ignore NA entries (which may arise from subsampling a knn-graph)
  if (all(is.na(knn.set))){ #if all values of a neighbourhood are NA
    return(NA)
  }
  else{
    freq.env <- table(batch[knn.set[!is.na(knn.set)]])/length(!is.na(knn.set))
    full.classes <- rep(0, length(class.freq$class))
    full.classes[ class.freq$class %in% names(freq.env)] <- freq.env
    exp.freqs <- class.freq$freq
    #compute chi-square test statistics
    sum((full.classes - exp.freqs)^2/exp.freqs)
  }
}

#which batch has the largest deviance (and is underrepresented)
max_deviance_batch <- function(knn.set, class.freq, batch) {
  #knn.set: indices of nearest neighbours
  #empirical frequencies in nn-environment (sample 1)
  if (all(is.na(knn.set))){ #if all values of a neighbourhood are NA
    return(NA)
  }
  else{
  freq.env <- table(batch[knn.set[!is.na(knn.set)]])/length(!is.na(knn.set))
  full.classes <- rep(0, length(class.freq$class))
  full.classes[ class.freq$class %in% names(freq.env)] <- freq.env
  exp.freqs <- class.freq$freq
  #compute chi-square test statistics
  allScores <- (full.classes - exp.freqs)/exp.freqs
  batch[which(allScores == min(allScores))]
  }
}


#the core function of kBET
chi_batch_test <- function(knn.set, class.freq, batch, df) {
  #knn.set: indices of nearest neighbours
  #empirical frequencies in nn-environment (sample 1)
  if (all(is.na(knn.set))){ #if all values of a neighbourhood are NA
    return(NA)
  }
  else{
  freq.env <- table(batch[knn.set[!is.na(knn.set)]])
  full.classes <- rep(0, length(class.freq$class))
  full.classes[ class.freq$class %in% names(freq.env)] <- freq.env
  exp.freqs <- class.freq$freq*length(knn.set)
  #compute chi-square test statistics
  chi.sq.value <- sum((full.classes - exp.freqs)^2/exp.freqs)
  result <- 1 - pchisq(chi.sq.value, df) #p-value for the result
  if (is.na(result)) { #I actually would like to now when 'NA' arises.
    return(NA)
  } else {
    result
  }
  }
}

lrt_approximation <- function(knn.set, class.freq, batch, df) {
  #knn.set: indices of nearest neighbours
  #empirical frequencies in nn-environment (sample 1)
  if (all(is.na(knn.set))){ #if all values of a neighbourhood are NA
    return(NA)
  }
  else{
  obs.env <- table(batch[knn.set[!is.na(knn.set)]]) #observed realisations of each category
  freq.env <- obs.env/sum(obs.env) #observed 'probabilities'
  full.classes <- rep(0, length(class.freq$class))
  obs.classes <- class.freq$class %in% names(freq.env)
  #for stability issues (to avoid the secret division by 0): introduce
  #another alternative model where the observed probability
  #is either the empirical frequency or 1/(sample size) at minimum
  if (length(full.classes) > sum(obs.classes)) {
    dummy.count <- length(full.classes) - sum(obs.classes)
    full.classes[obs.classes] <- obs.env / (sum(obs.env) + dummy.count)
    pmin <- 1/(sum(obs.env) + dummy.count)
    full.classes[!obs.classes] <- pmin
  } else {
    full.classes[ obs.classes] <- freq.env
  }
  exp.freqs <- class.freq$freq #expected 'probabilities'
  #compute likelihood ratio of null and alternative hypothesis,
  #test statistics converges to chi-square distribution
  full.obs <- rep(0, length(class.freq$class))
  full.obs[obs.classes] <- obs.env

  lrt.value <- -2*sum(full.obs * log(exp.freqs/full.classes))

  result <- 1 - pchisq(lrt.value, df) #p-value for the result
  if (is.na(result)) { #I actually would like to now when 'NA' arises.
    return(NA)
  } else {
    result
  }
  }
}
#truncated normal distribution distribution function
ptnorm <- function(x,mu,sd, a=0, b=1, alpha = 0.05, verbose = FALSE){
  #this is the cumulative density of the truncated normal distribution
  #x ~ N(mu, sd^2), but we condition on a <= x <= b
  if (a > b) {
    warning("Lower and upper bound are interchanged.")
    tmp <- a
    a <- b
    b <- tmp
  }

  if (sd <= 0 || is.na(sd)) {
    if (verbose) {
      warning("Standard deviation must be positive.")
    }
    if (alpha <= 0) {
      stop("False positive rate alpha must be positive.")
    }
    sd <- alpha
  }
  if (x < a || x > b) {
    warning("x out of bounds.")
    cdf <- as.numeric(x > a)
  } else {
    alp <- pnorm((a - mu) / sd)
    bet <- pnorm((b - mu) / sd)
    zet <- pnorm((x - mu) / sd)
    cdf <- (zet - alp) / (bet - alp)
  }
  cdf
}
#wrapper for the multinomial exact test function
multiNom <- function(x, y, z) {
  z.f <- factor(z)
  tmp <- multinomial.test(as.numeric(table(z.f[x])), y)
  tmp$p.value
}

#significance test for pcRegression (two levels)
correlate.fun_two <- function(rot.data, batch, batch.levels) {
  #rot.data: some vector (numeric entries)
  #batch: some vector (categoric entries)
  a <- lm(rot.data ~ batch)
  result <- numeric(2)
  result[1] <- summary(a)$r.squared #coefficient of determination
  result[2] <- summary(a)$coefficients[2,4] #p-value (significance level)
  t.test.result <- t.test(rot.data[batch == batch.levels[1]],
                          rot.data[batch == batch.levels[2]], paired = FALSE)
  result[3] <- t.test.result$p.value
  result
}

#significance test for pcRegression (more than two levels)
correlate.fun_gen <- function(rot.data, batch){
  #rot.data: some vector (numeric covariate)
  #batch: some vector (categoric covariate)
  a <- lm(rot.data ~ batch)
  result <- numeric(2)
  result[1] <- summary(a)$r.squared #coefficient of determination
  F.test.result <- aov(rot.data ~ batch)
  F.test.summary <- summary(F.test.result)

  result[2] <- summary(a)$coefficients[2,4] #p-value (significance level)
  result[3] <- F.test.summary[[1]]$'Pr(>F)'[1] #p-value of the one-way anova test

  result
}

