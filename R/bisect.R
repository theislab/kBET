#' bisect - a generic bisection function
#' @description provides recursive bisection algorithm for
#'    an arbitrary function
#'    It evaluates the function \code{foo} at the bounds and
#'    replaces one of the boundaries until a maximum is found
#'    or the interval becomes too small
#' @param foo a function mapping a one-dim argument to one-dim value
#'
#' @param bounds a vector of length 2 with real valued numbers
#'     (i.e. two arguments of \code{foo})
#' @param known tells for which of the arguments a value is known
#'     (defaults to NULL)
#' @param ... additional parameters for \code{foo}
#' @param tolx break condition for argument (defaults to 10)
#' @param toly break condition for value (defaults to 0.01)
#' @return A range of bounds where \code{foo} is maximal.
#' @importFrom stats dist
#'
#' @examples
#'     get_maximum <- bisect(function(x){-(x-2)^2}, c(-5,50))
#'
#' @export
bisect <- function(foo, bounds, known=NULL, ..., tolx = 5, toly = 0.01){
  if (is.null(known)) {
    evalFoo <- sapply(bounds, foo, ...)
    if (diff(unlist(evalFoo)) < -toly && diff(bounds) > tolx) {
      #i.e. the second bound has a smaller argument
      #than the first (lower) bound
      bounds[2] <- round(sum(bounds)/2,0)
      known <- c(unlist(evalFoo)[1],0)

      bisect(foo, bounds, known, ..., tolx = tolx, toly = toly)
    } else if (diff(unlist(evalFoo)) > toly && diff(bounds) > tolx) {
      #i.e. the first bound has a smaller argument
      #than the second (upper) bound
      bounds[1] <- round(sum(bounds)/2,0)
      known <- c(0,unlist(evalFoo)[2])

      bisect(foo, bounds,known,..., tolx = tolx, toly = toly)
    } else if (dist(unlist(evalFoo)) < toly) {
      #check the value in between upper and lower bound
      center <- sapply(round(sum(bounds)/2,0), foo, ...)
      dist_center <- sapply(unlist(evalFoo),
                            function(x,y){dist(c(x,y))},
                            center)
      if (max(abs(dist_center)) < toly || dist(bounds) < tolx) {
        #return interval and stop recursion
        bounds
      } else if (dist_center[1] > toly) {
        bounds[2] <- round(sum(bounds)/2,0)
        known <- c(unlist(evalFoo)[1],0)

        bisect(foo, bounds,known, ..., tolx = tolx, toly = toly)
      } else if (dist_center[2] > toly) {
        bounds[1] <- round(sum(bounds)/2,0)
        known <- c(0,unlist(evalFoo)[2])

        bisect(foo, bounds,known, ..., tolx = tolx, toly = toly)
      }
    }
  } else {
    new.eval <- which(known == 0)
    old.eval <- which(known != 0)
    evalFoo <- sapply(bounds[new.eval], foo, ...)
    result <- numeric(length(known))
    result[new.eval] <- unlist(evalFoo)
    result[old.eval] <- known[old.eval]
    if (diff(result) < -toly && diff(bounds) > tolx) {
      #i.e. the second bound has a smaller argument
      #than the first (lower) bound
      bounds[2] <- round(sum(bounds)/2,0)
      known <- c(unlist(evalFoo)[1],0)

      bisect(foo, bounds, known,..., tolx = tolx, toly = toly)
    } else if (diff(result) > toly && diff(bounds) > tolx) {
      #i.e. the first bound has a smaller argument
      #than the second (upper) bound
      bounds[1] <- round(sum(bounds)/2,0)
      known <- c(0,unlist(evalFoo)[1])

      bisect(foo, bounds, known,..., tolx = tolx, toly = toly)
    } else if (dist(result) < toly || dist(bounds) < tolx) {
      #return interval and stop recursion
      bounds
    }
  }
}
