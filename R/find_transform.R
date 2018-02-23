find_transform <- function(source_image,
                           target_image,
                           max_iter = 50,
                           inner_iter = 10,
                           max_search = 5,
                           window_size = 15,
                           resolution = 12,
                           scheme = 0,
                           tol = 1.e-5,
                           method = c("interpolation", "tps")) {

  patch_size <- window_size
  search_size <- max_search
  stride <- min(round(dim(source_image) / resolution)) #15
  sigma <- resolution / 4
  kernel <- gaussian_kernel(sigma = sigma, size = patch_size / 2)

  LL <- list()

  for (i in seq(inner_iter)) {

    LL[[i]] <- obtain_candidates_similarities(image = target_image,
                                              templates = array(source_image, dim = c(dim(target_image), 1)),
                                              patch_size = patch_size,
                                              search_size = search_size,
                                              stride = stride,
                                              max_iter = max_iter,
                                              max_random_neighbours = 5,
                                              crs = FALSE,
                                              early_stopping = list(tol = tol),
                                              scheme = scheme,
                                              ncores = parallel::detectCores() - 1)

  }

  my_voxels <- which(LL[[1]]$voxel_lookup_table >= 0)
  idx <- LL[[1]]$voxel_lookup_table[my_voxels] + 1

  R <- sapply(LL, function(x) x$similarities)
  cand <- sapply(LL, function(x) x$candidates)
  best_R <- apply(R, 1, which.min)
  my_candidates <- cand[LL[[1]]$actual_voxels * (best_R - 1) + seq_along(best_R)]
  my_similarities <- R[LL[[1]]$actual_voxels * (best_R - 1) + seq_along(best_R)]

  target_coords <- arrayInd(my_candidates[idx] + 1, .dim = dim(target_image)) - 1
  source_coords <- arrayInd(my_voxels, .dim = dim(target_image)) - 1


  transform <- compute_deformation_field(source_coords = source_coords,
                                         target_coords = target_coords,
                                         target_dims = dim(target_image),
                                         method = method[1],
                                         sigma = sigma,
                                         kernel = kernel)


  mapped_image <- deform_volume(V = source_image,
                                Dx = transform$Dx,
                                Dy = transform$Dy,
                                Dz = transform$Dz,
                                target_dims = dim(target_image),
                                method = 3)

  return(list(transform = transform,
              image = mapped_image,
              actual_voxels = LL[[1]]$actual_voxels,
              voxel_lookup_table = LL[[1]]$voxel_lookup_table,
              candidates = my_candidates,
              patch_neighbours = LL[[1]]$patch_neighbours,
              similarities = my_similarities))

}

compute_deformation_field <- function(source_coords,
                                      target_coords,
                                      target_dims,
                                      method = c("interpolation", "tps", "mixture"),
                                      sigma = 1,
                                      kernel_size = 0) {

  displ <- target_coords - source_coords

  x_unique <- unique(source_coords[, 1])
  y_unique <- unique(source_coords[, 2])
  z_unique <- unique(source_coords[, 3])
  stride <- x_unique[2] - x_unique[1]

  xout <- seq(min(x_unique), max(x_unique) + stride)
  yout <- seq(min(y_unique), max(y_unique) + stride)
  zout <- seq(min(z_unique), max(z_unique) + stride)

  SM <- scale_matrix(sx = 1 / stride,
                     sy = 1 / stride,
                     sz = 1 / stride)

  expanded_coords <- rbind(expand.grid(x = 0, y = y_unique, z = z_unique),
                           expand.grid(x = target_dims[1], y = y_unique, z = z_unique),
                           expand.grid(x = x_unique, y = 0, z = z_unique),
                           expand.grid(x = x_unique, y = target_dims[2], z = z_unique),
                           expand.grid(x = x_unique, y = y_unique, z = 0),
                           expand.grid(x = x_unique, y = y_unique, z = target_dims[3]))
  expanded_displ <- rep(0, times = nrow(expanded_coords))
  extended_df <- data.frame(x = expanded_coords[, 1],
                            y = expanded_coords[, 2],
                            z = expanded_coords[, 3],
                            D = expanded_displ)

  if (method[1] == "interpolation") {

    displ_array_x <- array(0, dim = c(length(x_unique) + 1, length(y_unique) + 1, length(z_unique) + 1))
    displ_array_x[seq_along(x_unique), seq_along(y_unique), seq_along(z_unique)] <- displ[, 1]

    Dx1 <- transform_volume(V = displ_array_x,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dx2 <- Dx1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    if (kernel_size > 0) {

      Dx_norm <- iterative_blur(image = Dx2, sigma = sigma, kernel_size = kernel_size, ncores = parallel::detectCores() - 1)

    } else Dx_norm <- Dx2

    displ_array_y <- array(0, dim = c(length(x_unique) + 1, length(y_unique) + 1, length(z_unique) + 1))
    displ_array_y[seq_along(x_unique), seq_along(y_unique), seq_along(z_unique)] <- displ[, 2]

    Dy1 <- transform_volume(V = displ_array_y,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dy2 <- Dy1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    if (kernel_size > 0) {

      Dy_norm <- iterative_blur(image = Dy2, sigma = sigma, kernel_size = kernel_size, ncores = parallel::detectCores() - 1)

    } else Dy_norm <- Dy2

    displ_array_z <- array(0, dim = c(length(x_unique) + 1, length(y_unique) + 1, length(z_unique) + 1))
    displ_array_z[seq_along(x_unique), seq_along(y_unique), seq_along(z_unique)] <- displ[, 3]
    Dz1 <- transform_volume(V = displ_array_z,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dz2 <- Dz1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    if (kernel_size > 0) {

      Dz_norm <- iterative_blur(image = Dz2, sigma = sigma, kernel_size = kernel_size, ncores = parallel::detectCores() - 1)

    } else Dz_norm <- Dz2

    return(list(Dx = Dx_norm, Dy = Dy_norm, Dz = Dz_norm))

  }

  if (method[1] == "tps") {

    require(mgcv)

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 1])
    # edf <- rbind(df, extended_df)
    b <- gam(D ~ te(x, y, z), data = df)
    # b <- gam(D ~ x * y * z, data = df)
    # displ_array_x <- array(predict(b, newdata = df), dim = c(length(x_unique), length(y_unique), length(z_unique)))
    displ_array_x <- array(b$fitted.values, dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dx1 <- transform_volume(V = displ_array_x,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dx2 <- Dx1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 2])
    # edf <- rbind(df, extended_df)
    b <- gam(D ~ te(x, y, z), data = df)
    # b <- gam(D ~ x * y * z, data = df)
    # displ_array_y <- array(unname(predict(b, newdata = df)), dim = c(length(x_unique), length(y_unique), length(z_unique)))
    displ_array_y <- array(b$fitted.values, dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dy1 <- transform_volume(V = displ_array_y,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dy2 <- Dy1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 3])
    # edf <- rbind(df, extended_df)
    b <- gam(D ~ te(x, y, z), data = df)
    # b <- gam(D ~ x * y * z, data = df)
    # displ_array_z <- array(unname(predict(b, newdata = df)), dim = c(length(x_unique), length(y_unique), length(z_unique)))
    displ_array_z <- array(b$fitted.values, dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dz1 <- transform_volume(V = displ_array_z,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dz2 <- Dz1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    return(list(Dx = Dx2, Dy = Dy2, Dz = Dz2))

  }

  if (method[1] == "poly") {

    degree <- 5

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 1])
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = TRUE), data = df)
    displ_array_x <- array(b$fitted.values, dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dx1 <- transform_volume(V = displ_array_x,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dx2 <- Dx1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 2])
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = TRUE), data = df)
    displ_array_y <- array(b$fitted.values, dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dy1 <- transform_volume(V = displ_array_y,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dy2 <- Dy1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 3])
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = TRUE), data = df)
    displ_array_z <- array(b$fitted.values, dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dz1 <- transform_volume(V = displ_array_z,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dz2 <- Dz1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    return(list(Dx = Dx2, Dy = Dy2, Dz = Dz2))

  }

  if (method[1] == "extended_poly") {

    degree <- 5

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 1])
    edf <- rbind(df, extended_df)
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = TRUE), data = edf)
    displ_array_x <- array(predict(b, newdata = df), dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dx1 <- transform_volume(V = displ_array_x,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dx2 <- Dx1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 2])
    edf <- rbind(df, extended_df)
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = TRUE), data = edf)
    displ_array_y <- array(predict(b, newdata = df), dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dy1 <- transform_volume(V = displ_array_y,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dy2 <- Dy1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 3])
    edf <- rbind(df, extended_df)
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = TRUE), data = edf)
    displ_array_z <- array(predict(b, newdata = df), dim = c(length(x_unique), length(y_unique), length(z_unique)))

    Dz1 <- transform_volume(V = displ_array_z,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dz2 <- Dz1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    return(list(Dx = Dx2, Dy = Dy2, Dz = Dz2))

  }

  if (method[1] == "filterpoly") {

    degree <- 5
    raw <- FALSE
    rs <- rowSums(abs(displ))
    q <- quantile(rs, probs = c(0.25, 0.75))
    index <- which(rs >= q[2] | rs <= q[1])

    new_source <- source_coords[index, ]
    new_target <- target_coords[index, ]
    new_displ <- displ[index, ]


    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 1])
    new_df <- data.frame(x = new_source[, 1], y = new_source[, 2], z = new_source[, 3], D = new_displ[, 1])
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = raw), data = new_df)
    L <- sapply(apply(source_coords, 2, unique), length)
    displ_array_x <- array(predict(b, newdata = df), dim = L)

    Dx1 <- transform_volume(V = displ_array_x,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dx2 <- Dx1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 2])
    new_df <- data.frame(x = new_source[, 1], y = new_source[, 2], z = new_source[, 3], D = new_displ[, 2])
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = raw), data = new_df)
    displ_array_y <- array(predict(b, newdata = df), dim = L)

    Dy1 <- transform_volume(V = displ_array_y,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dy2 <- Dy1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 3])
    new_df <- data.frame(x = new_source[, 1], y = new_source[, 2], z = new_source[, 3], D = new_displ[, 3])
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = raw), data = new_df)
    displ_array_z <- array(predict(b, newdata = df), dim = L)

    Dz1 <- transform_volume(V = displ_array_z,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dz2 <- Dz1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    return(list(Dx = Dx2, Dy = Dy2, Dz = Dz2))



  }

  if (method[1] == "extfiltpoly") {

    degree <- 5
    raw <- FALSE
    rs <- rowSums(abs(displ))
    probs <- 0.35
    q <- quantile(rs, probs = c(probs, 1 - probs))
    index <- which(rs >= q[2] | rs <= q[1])

    new_source <- source_coords[index, ]
    new_target <- target_coords[index, ]
    new_displ <- displ[index, ]


    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 1])
    new_df <- data.frame(x = new_source[, 1], y = new_source[, 2], z = new_source[, 3], D = new_displ[, 1])
    new_df <- rbind(new_df, extended_df)
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = raw), data = new_df)
    L <- sapply(apply(source_coords, 2, unique), length)
    displ_array_x <- array(predict(b, newdata = df), dim = L)

    Dx1 <- transform_volume(V = displ_array_x,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dx2 <- Dx1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 2])
    new_df <- data.frame(x = new_source[, 1], y = new_source[, 2], z = new_source[, 3], D = new_displ[, 2])
    new_df <- rbind(new_df, extended_df)
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = raw), data = new_df)
    displ_array_y <- array(predict(b, newdata = df), dim = L)

    Dy1 <- transform_volume(V = displ_array_y,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dy2 <- Dy1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    df <- data.frame(x = source_coords[, 1], y = source_coords[, 2], z = source_coords[, 3], D = displ[, 3])
    new_df <- data.frame(x = new_source[, 1], y = new_source[, 2], z = new_source[, 3], D = new_displ[, 3])
    new_df <- rbind(new_df, extended_df)
    b <- lm(D ~ poly(x, y, z, degree = degree, raw = raw), data = new_df)
    displ_array_z <- array(predict(b, newdata = df), dim = L)

    Dz1 <- transform_volume(V = displ_array_z,
                            M = SM,
                            target_dims = c(length(xout), length(yout), length(zout)),
                            method = 3)

    Dz2 <- Dz1[seq(target_dims[1]), seq(target_dims[2]), seq(target_dims[3])]

    return(list(Dx = Dx2, Dy = Dy2, Dz = Dz2))

  }

  # Default: no transformation
  return(list(Dx = array(0, dim = target_dims), Dy = array(0, dim = target_dims), Dz = array(0, dim = target_dims)))

}
