#' Transform the Graylevels of an Image by Predefined Intervals
#'
#' @param input                   (3D image) Image to be transformed
#' @param input_intervals         (Numeric vector) Ordered values which are considered 
#' graylevel intervals in the input image. See Details.
#' @param reference_intervals     (Numeric vector) Ordered values to which map the input intervals.
#'
#' @details Suppose \code{input_intervals = c(0, 10, 25, 40)} (0 and 40 are the minimum and maximum of the \code{input} image, respectively) and \code{reference_intervals = c(5, 70, 90, 100)}. This transformation linearly maps, consecutively, interval [0, 10] to [5, 70], then [10, 25] to [70, 90], and so.
#'
#' @return The transformed image.
#'
map_intervals <- function(input, input_intervals, reference_intervals, verbose = FALSE) {
  
  # Number of intervals
  n_intervals <- length(input_intervals) - 1
  
  # For each interval, select the gray-levels that lie withtin its extreme values and map their values to 
  # the corresponding reference interval, with a linear interpolation
  res <- 0 * input
  
  for (i in 1:n_intervals) {
    
    input_lo <- input_intervals[i]
    input_hi <- input_intervals[i + 1]
    
    ref_lo <- reference_intervals[i]
    ref_hi <- reference_intervals[i + 1]
    
    if (verbose)
      cat("[", input_lo, ",", input_hi, "] -> [", ref_lo, ",", ref_hi, "]\n")
    
    idx <- which(input >= input_lo & input <= input_hi)
    res[idx] <- (input[idx] - input_lo) / (input_hi - input_lo) * (ref_hi - ref_lo) + ref_lo
    
  }
  
  return(res)
  
}