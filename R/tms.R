#' Robust Estimation of Average Gray-Level of GM, WM and CSF
#'
#' @param img    (3D image) Image to segment
#' @param k      (Integer) Number of classes to segment. Currently, only 3.
#'
#' @return A vector with 3 components, indicating the estimated mean gray-level for each tissue.
#'
tms <- function(img, k = 3) {
  
  alfa <- 0.05 # 5% data of smaller gradient magnitude
  
  # Smooth the image
  ima <- cvolfilter2d(ima = as.array(img), dims = dim(img), v = 1)
  
  # Obtain gradient of the image
  vima <- cgradientmodule(ima = ima, dims = dim(ima))
  
  # Select points of the image with little gradient (non-edges) and perform 
  # recursive k-means to compute the 3 classes
  ind <- which(ima > 0)
  m <- sd(vima[ind])
  
  mascara <- rkmeans(ima, mask = (vima <= 2 * m))
  mascara <- cTruncado(mascara$mask, dims = dim(mascara$mask))
  
  # For each class, compute the mean of its gray-level using least trimmed squares approximation
  loc <- rep(0, times = k)
  
  for (i in 1:k) {
    
    ind <- which(mascara == i)
    p <- ima[ind]  
    v <- vima[ind]
    
    h <- histogram(v)
    
    t2 <- sum(h) * alfa
    
    h2 <- cumsum(h)
    
    n <- which(h2 > t2)[1]
    
    ind <- which(v < n)
    p <- p[ind]
    
    loc[i] <- lts(p)
    
  }
  
  c <- sort(loc)

  return(c)
  
}





