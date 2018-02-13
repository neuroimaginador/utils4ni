#' Normalize the Image to a Reference Histogram
#'
#' @param img           (3D array) Image to normalize.
#' @param tms_params    (Numeric vector) Vector with 3 components, indicating the mean values of the 3 tissues. If not provided, this function computes them using \code{\link{tms}}.
#' @param reference_intervals   (Numeric vector) Segments of the histogram. Default: c(0, 50, 150, 250, 255). Must have 5 elements.
#'
#' @return The transformed image
#' 
intensity_normalization_Otsu <- function(img, mask = img > 0, reference_intervals = c(0, 50, 150, 250, 255)) {
  
  require(mritc)
  
  tms_params <- get_contrast_estimation(img, mask, method = "otsu")
  
  res <- map_intervals(input = as.array(img), 
                       input_intervals = c(min(img), tms_params, max(img)),
                       reference_intervals = reference_intervals)
  
  return(list(params = tms_params, image = res))
  
}