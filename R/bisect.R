
#create bisection function
#it evaluates the function foo at the bounds (a vector with exactly 2 entries) and
#replaces one of the boundaries until a maximum is found or the interval becomes too small
#' @export
bisect <- function(foo, bounds, known=NULL,..., tolx = 10, toly = 0.01){
  if(is.null(known)){
    evalFoo <- sapply(bounds, foo, ...)
    if (diff(unlist(evalFoo))< -toly& diff(bounds)> tolx){ #i.e. the second bound has a smaller argument than the first (lower) bound
      bounds[2] <- round(sum(bounds)/2,0)
      known <- c(unlist(evalFoo)[1],0)
      return(bisect(foo, bounds,known,..., tolx = tolx, toly = toly))
    }
    if (diff(unlist(evalFoo))>toly& diff(bounds)> tolx){ #i.e. the first bound has a smaller argument than the second (upper) bound
      bounds[1] <- round(sum(bounds)/2,0)
      known <- c(0,unlist(evalFoo)[2])
      return(bisect(foo, bounds,known,..., tolx = tolx, toly = toly))
    }
    if(dist(unlist(evalFoo))<toly | dist(bounds) < tolx){
      #return interval and stop recursion
      return(bounds)
    }
  }else{
    new.eval <- which(known==0)
    old.eval <- which(known!=0)
    evalFoo <- sapply(bounds[new.eval], foo, ...)
    result <- numeric(length(known))
    result[new.eval] <- unlist(evalFoo)
    result[old.eval] <- known[old.eval]
    if (diff(result)< -toly & diff(bounds)> tolx){ #i.e. the second bound has a smaller argument than the first (lower) bound
      bounds[2] <- round(sum(bounds)/2,0)
      known <- c(unlist(evalFoo)[1],0)
      return(bisect(foo, bounds, known,..., tolx = tolx, toly = toly))
    }
    if (diff(result)>toly & diff(bounds)> tolx){ #i.e. the first bound has a smaller argument than the second (upper) bound
      bounds[1] <- round(sum(bounds)/2,0)
      known <- c(0,unlist(evalFoo)[1])
      return(bisect(foo, bounds, known,..., tolx = tolx, toly = toly))
    }
    if(dist(result)<toly | dist(bounds) < tolx){
      #return interval and stop recursion
      return(bounds)
    }
  }

}
