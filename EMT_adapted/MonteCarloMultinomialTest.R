MonteCarloMultinomialTest <-function (observed, prob, size, groups, numEvents, ntrial, atOnce, verbose) 
{
  IDofObs <- paste(observed, sep = "", collapse = "")
  sumLowFreq = 0
  loops = floor(ntrial/atOnce)
  rest = ntrial%%atOnce
  if (loops > 0) {
    for (i in 1:loops) {
      res <- rmultinom(n = atOnce, size = size, prob = prob)
      vec <- apply(res, 2, function(x) paste(x, sep = "", 
                                             collapse = ""))
      frequencyTable <- table(vec)
      if (sum(rownames(frequencyTable) == IDofObs) == 0) {
        freqObs = 0
      }
      else {
        freqObs <- frequencyTable[rownames(frequencyTable) == 
                                    IDofObs][[1]]
      }
      sumLowFreq <- sumLowFreq + sum(frequencyTable[as.vector(frequencyTable) <= 
                                                      freqObs])
      cat(" Number of withdrawals accomplished: ", prettyNum(i * 
                                                               atOnce, scientific = FALSE, big.mark = "."), 
          "\n")
      flush.console()
    }
  }
  if (rest > 0) {
    res <- rmultinom(n = rest, size = size, prob = prob)
    vec <- apply(res, 2, function(x) paste(x, sep = "", collapse = ""))
    frequencyTable <- table(vec)
    if (sum(rownames(frequencyTable) == IDofObs) == 0) {
      freqObs = 0
    }
    else {
      freqObs <- frequencyTable[rownames(frequencyTable) == 
                                  IDofObs][[1]]
    }
    sumLowFreq <- sumLowFreq + sum(frequencyTable[as.vector(frequencyTable) <= 
                                                    freqObs])
  }
  p.value = (sumLowFreq + 1)/(ntrial + 1)
  if (verbose){
  head <- paste("\n Monte Carlo Multinomial Test, distance measure: f\n\n")
  tab <- as.data.frame(cbind(numEvents, round(freqObs/ntrial, 
                                              digits = 4), round(p.value, digits = 4)))
  colnames(tab) <- c("   Events", "   fObs", "   p.value")
  cat(head)
  print(tab, row.names = FALSE)
  }
  invisible(list(id = "Monte Carlo Multinomial Test", size = size, 
                 groups = groups, stat = "lowF", allProb = sort(as.vector(frequencyTable), 
                                                                decreasing = TRUE)/ntrial, ntrial = ntrial, p.value = round(p.value, 
                                                                                                                            digits = 4)))
}