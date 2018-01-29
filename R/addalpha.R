#' addalpha - add transparency info to colorset for plotting
#' @description \code{addalpha} adds transparency
#'     information to a vector in R color format
#' @param colors a hexadecimal string as commonly used in R
#'
#' @param alpha
#'     a transparency factor:
#'     0 - completely invisible
#'     1 - fully opaque
#'
#' @return a hexadecimal string of length \code{colors}
#' @examples
#'     library(RColorBrewer)
#'     set2 <- brewer.pal(5, 'Set2')
#'     set2.light <- addalpha(set2, alpha=0.5)
#'
#' @importFrom grDevices col2rgb rgb
#' @export
addalpha <- function(colors, alpha=1.0) {
    r <- col2rgb(colors, alpha=TRUE)
    # Apply alpha
    r[4,] <- alpha*255
    r <- r/255.0
    return(rgb(r[1,], r[2,], r[3,], r[4,]))
}
