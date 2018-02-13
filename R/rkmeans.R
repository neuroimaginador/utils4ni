#' Recursive 3-Means for Image Segmentation
#'
#' @param imag     (3D array o coerceable to it) Image to segment into 3 classes.
#' @param mask     (3D array) Mask with the region of interest inside the image
#'
#' @return A list with two elements: \code{means}, the mean gray level of the 3 classes; and \code{mask}, a binary mask of the segmented classes.
#'
rkmeans <- function(imag, mask) {
  
  # Image histogram inside the mask
  histo <- histogram(imag[mask])
  
  # Mean value of the image inside the mask, according to the histogram
  ind <- which(histo > 0)
  mu <- sum(histo[ind] * ind) / sum(histo[ind])
  
  # Find best classification into 2 classes
  mu <- c(mu - 1, mu + 1)
  means <- findmeans(histo, mu)

  # Find best classification into 3 classes
  m <- c(means$mu[1] + 1, means$mu[1] - 1, means$mu[2])
  bestmeans1 <- findmeans(histo, m)
  
  m <- c(means$mu[2] + 1, means$mu[2] - 1, means$mu[1])
  bestmeans2 <- findmeans(histo, m)
  
  bestmu <- bestmeans1$mu
  if (bestmeans2$error < bestmeans1$error) {
    
    bestmu <- bestmeans2$mu
    
  }
  
  # Create mask
  mask <- cCreateMask(imag, dims = dim(imag), c = 3, medias = sort(bestmu))
  
  return(list(means = bestmu, mask = mask))
  
}
