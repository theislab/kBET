#from EMT package (, function is needed to run ExactMultinomialTest.R)

findVectors <- function (groups, size)
{
  if (groups == 1) {
    mat = size
  }
  else {
    mat <- matrix(rep(0, groups - 1), nrow = 1)
    for (i in seq_len(size)) {
      mat <- rbind(mat, findVectors(groups - 1, i))
    }
    mat <- cbind(mat, size - rowSums(mat))
  }
  invisible(mat)
}
