#' Histogram Computation
#'
#' @param data    (a vector)  Data to extract estimate the histogram from.
#'
#' @return A vector with the densities from building the histogram of the \code{data}.
#' 
histogram <- function(data) {

  M <- ceiling(max(data)) + 1
  
  h <- hist(data, breaks = M, plot = FALSE)$density

  return(h)

}
