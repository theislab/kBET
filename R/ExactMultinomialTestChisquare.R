#from EMT package (solemnly adapted to stay quiet)

ExactMultinomialTestChisquare <- function (observed, prob, size, groups, numEvents, verbose)
{
  expectedFreq = size * prob
  chi2Obs = chisqStat(observed, expectedFreq)
  eventMat <- findVectors(groups, size)
  if (nrow(eventMat) != numEvents) 
    stop("Wrong number of events calculated. \n This is probably a bug.")
  eventProb <- apply(eventMat, 1, function(x) dmultinom(x, 
                                                        size = size, prob = prob))
  eventChi2 <- apply(eventMat, 1, function(x) chisqStat(x, 
                                                        expectedFreq))
  eventPandChi2 <- cbind(eventProb, eventChi2)
  if (round(sum(eventProb), digits = 2) != 1) 
    stop("Wrong values for probabilities. \n This is probably a bug.")
  p.value <- sum((eventPandChi2[eventPandChi2[, 2] >= chi2Obs, 
                                ])[, 1])
  if (verbose){
  head <- paste("\n Exact Multinomial Test, distance measure: chisquare\n\n")
  tab <- as.data.frame(cbind(numEvents, round(chi2Obs, digits = 4), 
                             round(p.value, digits = 4)))
  colnames(tab) <- c("   Events", "   chi2Obs", "   p.value")
  cat(head)
  print(tab, row.names = FALSE)
  }
  invisible(list(id = "Exact Multinomial Test, Chisquare", 
                 size = size, groups = groups, stat = "highChisq", allProb = sort(eventProb, 
                                                                                  decreasing = TRUE), ntrial = NULL, p.value = round(p.value, 
                                                                                                                                     digits = 4)))
}
