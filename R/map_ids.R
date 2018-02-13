map_ids_cpp <- function(image, remap_classes, invert = FALSE) {

  s <- remap_classes$source
  t <- remap_classes$target
  if (!invert) {

    Y <- map_ids_workhorse(x = image, source = s, target = t)


    remaining <- 0

    if (!is.null(remap_classes$remaining))
      remaining <- remap_classes$remaining

    Y2 <- map_extra_classes(x = image, source = s, remaining = remaining)
    Y2[Y > 0] <- 0
    Y <- Y + Y2

  } else {

    Y <- map_ids_workhorse(x = image, source = t, target = s)

  }

  return(Y)

}
