get_image_stats <- function(V) {

  v <- as.vector(V)

  stats <- list(
    min = min(v),
    max = max(v),
    mean = mean(v),
    sd = sd(v)
  )

  return(stats)

}
