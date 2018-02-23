iterative_blur <- function(A, sigma, w = 3) {

  n <- 3 * round(sigma)
  kernel <- gaussian_kernel(sigma = sigma, size = w)
  # kernel <- array(1 / (w ^ 3), dim = c(w, w, w))

  for (i in seq(n)) A <- regularize(A, kernel)

  return(A)

}
