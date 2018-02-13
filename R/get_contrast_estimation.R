#' @title Compute Parameters for Contrast Estimation
#' This function computes some parameters that may help in estimating the actual contrast in an image.
#'
#' @param img     (nifti object or 3D array) The brain image to be analyzed. Must be brain-extracted (skull stripped).
#'
#' @return A list with two elements: \code{params}, a vector with 3 components, which estimate the mean gray-level for each tissue type (CSF, GM and WM); and \code{contrast_coefficient}, a number indicating the contrast.
#'
get_contrast_estimation <- function(img, mask = img > 0, method = c("otsu", "tms")) {
  
  # Obtain a estimation of the values of the gray levels for CSF, GM and WM
  if (method == "otsu") {
    
    require(mritc)
    
    y <- img[mask > 0]
    
    m <- min(y)
    M <- max(y)
    y <- (y - m) / (M - m) * 512
    y <- round(y)
    params <- initOtsu(y, 2)$mu / 512 * (M - m) + m
    
  } else {
    
    params <- tms(img, k = 3)
    
  }
  
  return(list(params = params, contrast_coefficient = params[3] / params[2]))
  
}