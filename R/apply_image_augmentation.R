apply_image_augmentation <- function(V, M, type = c("continuous", "categorical")) {

  if (all(M == identity_matrix())) return(V)

  interp_method <- switch(tolower(type[1]),
                          "continuous" = 3,
                          "categorical" = 0)

  n_dims <- length(dim(V))

  tmpV <- V

  if (n_dims == 4) {

    num_volumes <- dim(V)[4]
    volume_dims <- dim(V)[1:3]

    for (i in seq(num_volumes)) {

      tmpV[, , , i] <- transform_volume(V = V[, , , i],
                                        M = M,
                                        target_dims = volume_dims,
                                        method = interp_method)

    }

  } else {

    tmpV <- transform_volume(V = V,
                             M = M,
                             target_dims = dim(V),
                             method = interp_method)

  }

  return(tmpV)

}
