#' Least Trimmed Squares for Robust Class Mean Computation
#'
#' @param y     (Numeric vector) Vector to compute its mean.
#' @param h     (Integer) Window size for the rolling mean computations
#'
#' @return The mean of the vector.
#'
lts <- function(y, h = NULL) {

  # We use a trick to compute the LTS fit of the mean
  # First, we sort the values
  y <- base::sort(y)

  # Compute the rolling means of given width
  if (is.null(h)) {
    
    h <- max(1, round(length(y) / 2))
    
  } 
  
  means <- zoo::rollmean(y, k = h)
  meansquares <- zoo::rollmean(y ^ 2, k = h)

  # Compute an approach to the standard deviation, using these values
  sds <- meansquares - means ^ 2

  # Select the point with the lowest standard deviation and return it
  idx <- which.min(sds)

  return(means[idx])

}

rollmean <- function(x, k) {
  
  n <- length(x)
  
  stopifnot(k <= n)
  
  ix <- floor((1 + k)/2):ceiling(n - k/2)
  xu <- unclass(x)
  y <- xu[k:n] - xu[c(1, seq_len(n - k))]
  y[1] <- sum(xu[1:k])
  rval <- cumsum(y)/k
  x[ix] <- rval
  zoo:::na.fi
  
  return(x)
  
}

# 
# lts <- function(y, h) {
#   
#   mean(y, trim = 0.02)
#   
# }