#' Denoise an Image
#'
#' @param img   (3D array) Image to be denoised
#'
#' @return The denoised image.
#'
adaptive_denoising <- function(img) {
  
  # Use C code to perform the adaptive nonlocal means denoising algorithm
  res <- adaptive_nonlocal_means_denoising(ima = as.array(img), dims = dim(img), v = 3, f = 1, h = 2)
  
  return(res)
  
}
