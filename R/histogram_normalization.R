#' Normalize the Image to a Reference Histogram
#'
#' @param img                   (3D array) Image to normalize.
#' @param reference_intervals   (Numeric vector) Segments of the histogram. Default: c(0, 50, 150, 250, 255). Must have 5 elements.
#'
#' @return The transformed image.
#'
histogram_normalization <- function(img, reference_intervals = c(0, 50, 150, 250, 255)) {
  
  # Compute gray-level of the 3 classes
  TMS_params <- tms(img, k = 3)
  
  # Map the image with this intervales to the reference one.
  map_intervals(input = as.array(img), 
                input_intervals = c(min(img), TMS_params, max(img)),
                reference_intervals = reference_intervals)
  
  
}
