#from EMT package (solemnly adapted to stay quiet)

ExactMultinomialTest <-function (observed, prob, size, groups, numEvents, verbose) 
{
  pObs = dmultinom(observed, size = size, prob)
  eventMat <- findVectors(groups, size)
  if (nrow(eventMat) != numEvents) 
    stop("Wrong number of events calculated. \n This is probably a bug.")
  eventProb <- apply(eventMat, 1, function(x) dmultinom(x, 
                                                        size = size, prob = prob))
  p.value = sum(eventProb[eventProb <= pObs])
  if (round(sum(eventProb), digits = 2) != 1) 
    stop("Wrong values for probabilities. \n This is probably a bug.")
  
  
  if (verbose){
    head <- paste("\n Exact Multinomial Test, distance measure: p\n\n")
    tab <- as.data.frame(cbind(numEvents, round(pObs, digits = 4), 
                               round(p.value, digits = 4)))
    colnames(tab) <- c("   Events", "   pObs", "   p.value")
    cat(head)
    print(tab, row.names = FALSE)
  }
  
  invisible(list(id = "Exact Multinomial Test", size = size, 
                 groups = groups, stat = "lowP", allProb = sort(eventProb, 
                                                                decreasing = TRUE), ntrial = NULL, p.value = round(p.value, 
                                                                                                                   digits = 4)))
}