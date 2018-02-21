map <- function(source, source_mask = NULL, target, target_mask = NULL, nbins = 64) {

  kernel <- gaussian_kernel(sigma = 4, size = 5)
  s <- regularize(source, kernel, ncores = 3)
  t <- regularize(target, kernel, ncores = 3)

  s <- s / max(s) * 255
  t <- t / max(t) * 255

  # s <- source
  if (!is.null(source_mask)) s <- s[source_mask > 0]
  s <- as.vector(s)

  # t <- target
  if (!is.null(target_mask)) t <- t[target_mask > 0]
  t <- as.vector(t)

  mx <- min(s)
  Mx <- max(s)

  my <- min(t)
  My <- max(t)

  x <- floor(nbins * (s - mx) / (Mx - mx))
  y <- floor(nbins * (t - my) / (My - my))

  x_centers <- (seq(nbins + 1) - 1) / nbins * (Mx - mx) + mx
  y_centers <- (seq(nbins + 1) - 1) / nbins * (My - my) + my
  M <- confusion_matrix(x, y)

  # v <- apply(M, 1, which.max)

  v <- apply(M, 1, function(row) {

    sum(seq(nbins + 1) * row) / sum(row)

  })

  new_centers <- approx(seq(nbins + 1), y = y_centers, xout = v)$y
  #
  # new_centers <- c(0, y_centers[v])
  new_x <- approx(x = x_centers, y = new_centers, xout = source / max(source) * 255)$y
  dim(new_x) <- dim(source)

  new_x <- new_x / max(new_x) * 255

  return(new_x)

}