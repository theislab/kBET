MonteCarloMultinomialTestChisquare <- function (observed, prob, size, groups, numEvents, ntrial, atOnce) 
{
  expectedFreq = size * prob
  chi2Obs = chisqStat(observed, expectedFreq)
  bigChis = 0
  loops = floor(ntrial/atOnce)
  rest = ntrial%%atOnce
  if (loops > 0) {
    for (i in 1:loops) {
      res <- rmultinom(n = atOnce, size = size, prob = prob)
      chi2all <- apply(res, 2, function(x) chisqStat(x, 
                                                     expectedFreq))
      bigChis <- bigChis + length(chi2all[chi2all >= chi2Obs])
      cat(" Number of withdrawals accomplished: ", prettyNum(i * 
                                                               atOnce, scientific = FALSE, big.mark = "."), 
          "\n")
      flush.console()
    }
  }
  if (rest > 0) {
    res <- rmultinom(n = rest, size = size, prob = prob)
    chi2all <- apply(res, 2, function(x) chisqStat(x, expectedFreq))
    bigChis <- bigChis + length(chi2all[chi2all >= chi2Obs])
  }
  if ((atOnce * loops + rest) != ntrial) 
    stop(" Number of withdrawals made incorrect.\n This is probably a bug.")
  p.value = (bigChis + 1)/(ntrial + 1)
  if (verbose){
  head <- paste("\n Monte Carlo Multinomial Test, distance measure: chisquare\n\n")
  tab <- as.data.frame(cbind(numEvents, round(chi2Obs, digits = 4), 
                             round(p.value, digits = 4)))
  colnames(tab) <- c("   Events", "   chi2Obs", "   p.value")
  cat(head)
  print(tab, row.names = FALSE)
  }
  invisible(list(id = "Monte Carlo Multinomial Test, Chisquare", 
                 size = size, groups = groups, stat = "highChi2", allProb = sort(as.vector(table(chi2all))/ntrial, 
                                                                                 decreasing = TRUE), ntrial = ntrial, p.value = round(p.value, 
                                                                                                                                      digits = 4)))
}