iterative_blur <- function(image, sigma, kernel_size = 3, ncores = parallel::detectCores() - 1) {

  n <- 3 * round(sigma)
  kernel <- gaussian_kernel(sigma = sigma, size = kernel_size)
  # kernel <- array(1 / (w ^ 3), dim = c(w, w, w))

  A <- image
  for (i in seq(n)) A <- regularize(A, kernel, ncores = ncores)

  return(A)

}
