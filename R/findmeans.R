#' Find Best Class Means for Gray-Level Segmentation
#'
#' @param h     (Numeric vector) Histogram of an image.
#' @param mu    (Numeric vector) Candidate class means.
#'
#' @return A list with 2 components: \code{mu}, a vector with the fitted means; and \code{error} the intra-class squared sum of distances.
#'
findmeans <- function(h, mu) {
  
  lh <- length(h)
  k <- length(mu)
  hc <- rep(0, times = lh)
  
  bool <- TRUE
  while (bool) {
    
    oldmu <- mu
    
    for (i in 1:lh) {
      
      c <- abs(i - mu)
      hc[i] <- which.min(c)
      
    }
    
    for (i in 1:k) {
      
      a <- which(hc == i)
      
      if (sum(h[a]) > 0)  mu[i] <- sum(a * h[a]) / sum(h[a])
      
    }
    
    bool <- any(abs(mu - oldmu) > 0)
    # print(mu)
  }
  
  error <- 0
  for (i in 1:k) {
    
    a <- which(hc == i)
    if (length(a) > 0) {
      
      error <- error + sum(abs(a - mu[i]) * h[a])
      
    }
    
  }

  return(list(error = error, mu = mu))
  
}
