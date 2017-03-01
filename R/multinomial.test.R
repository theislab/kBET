#from EMT package (function needed to run ExactMultinomialTest)

multinomial.test <- function (observed, prob, useChisq = FALSE, MonteCarlo = FALSE, ntrial = 1e+05, atOnce = 1e+06, verbose=FALSE)
{
  if (!is.vector(observed, mode = "numeric"))
    stop(" Observations have to be stored in a vector, e.g.  'observed <- c(5,2,1)'")
  if (!is.vector(prob, mode = "numeric"))
    stop(" Probabilities have to be stored in a vector, e.g.  'prob <- c(0.25, 0.5, 0.25)'")
  if (round(sum(prob), digits = 1) != 1)
    stop("Wrong input: sum of probabilities must not deviate from 1.")
  if (length(observed) != length(prob))
    stop(" Observations and probabilities must have same dimensions.")
  size <- sum(observed)
  groups <- length(observed)
  numEvents <- choose(size + groups - 1, groups - 1)
  res <- ExactMultinomialTest(observed, prob, size,
                              groups, numEvents, verbose)

  invisible(res)
}
