
#create bisection function
#it evaluates the function foo at the bounds (a vector with exactly 2 entries) and
#replaces one of the boundaries until a maximum is found or the interval becomes too small
#' @export
bisect <- function(foo, bounds,..., tolx = 10, toly = 0.01){
  evalFoo <- sapply(bounds, foo, ...)
  if (diff(unlist(evalFoo))< -toly){ #i.e. the second bound has a smaller argument than the first (lower) bound
    bounds[2] <- round(sum(bounds)/2,0)
    return(bisect(foo, bounds,..., tolx = tolx, toly = toly))
  }
  if (diff(unlist(evalFoo))>toly){ #i.e. the first bound has a smaller argument than the second (upper) bound
    bounds[1] <- round(sum(bounds)/2,0)
    return(bisect(foo, bounds,..., tolx = tolx, toly = toly))
  }
  if(dist(unlist(evalFoo))<toly || dist(bounds) < tolx){
    #return interval and stop recursion
    return(bounds)
  }
}
