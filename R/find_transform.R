find_transform <- function(source_image,
                           target_image,
                           max_iter = 50,
                           inner_iter = 10,
                           max_search = 5,
                           window_size = 15,
                           resolution = 12,
                           scheme = 0,
                           tol = 1.e-5) {

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

  displ <- target_coords - source_coords

  x_unique <- unique(source_coords[, 1])
  y_unique <- unique(source_coords[, 2])
  z_unique <- unique(source_coords[, 3])
  stride <- x_unique[2] - x_unique[1]

  xout <- seq(min(x_unique), max(x_unique))
  yout <- seq(min(y_unique), max(y_unique))
  zout <- seq(min(z_unique), max(z_unique))

  SM <- scale_matrix(sx = 1 / stride,
                     sy = 1 / stride,
                     sz = 1 / stride)

  displ_array_x <- array(0, dim = c(length(x_unique), length(y_unique), length(z_unique)))
  displ_array_x[] <- displ[, 1]

  Dx1 <- transform_volume(V = displ_array_x,
                          M = SM,
                          target_dims = c(length(xout), length(yout), length(zout)),
                          method = 3)

  Dx2 <- array(0, dim = dim(target_image))

  offsets <- dim(Dx2) - dim(Dx1)
  x_init <- round(offsets[1] / 2)
  if (x_init < 1) x_init <- 1
  x_end <- x_init + dim(Dx1)[1] - 1
  y_init <- round(offsets[2] / 2)
  if (y_init < 1) y_init <- 1
  y_end <- y_init + dim(Dx1)[2] - 1
  z_init <- round(offsets[3] / 2)
  if (z_init < 1) z_init <- 1
  z_end <- z_init + dim(Dx1)[3] - 1


  Dx2[x_init:x_end,
      y_init:y_end,
      z_init:z_end] <- Dx1

  Dx_norm <- regularize(image = Dx2, kernel = kernel, ncores = parallel::detectCores() - 1)

  displ_array_y <- array(displ[, 2], dim = c(length(x_unique), length(y_unique), length(z_unique)))

  Dy1 <- transform_volume(V = displ_array_y,
                          M = SM,
                          target_dims = c(length(xout), length(yout), length(zout)),
                          method = 3)

  Dy2 <- array(0, dim = dim(target_image))
  Dy2[x_init:x_end,
      y_init:y_end,
      z_init:z_end] <- Dy1

  Dy_norm <- regularize(image = Dy2, kernel = kernel, ncores = parallel::detectCores() - 1)

  displ_array_z <- array(displ[, 3], dim = c(length(x_unique), length(y_unique), length(z_unique)))
  Dz1 <- transform_volume(V = displ_array_z,
                          M = SM,
                          target_dims = c(length(xout), length(yout), length(zout)),
                          method = 3)

  Dz2 <- array(0, dim = dim(target_image))
  Dz2[x_init:x_end,
      y_init:y_end,
      z_init:z_end] <- Dz1

  Dz_norm <- regularize(image = Dz2, kernel = kernel, ncores = parallel::detectCores() - 1)

  # V1 <- deform_volume(V = source_image, Dx2, Dy2, Dz2, target_dims = dim(target_image), method = 3)
  # # ortho_plot(V1, text = "Dx")
  V1_norm <- deform_volume(V = source_image, Dx_norm, Dy_norm, Dz_norm, target_dims = dim(target_image), method = 3)
  # ortho_plot(V1_norm, text = "Dx_norm")

  return(list(transform = list(Dx = Dx_norm,
                               Dy = Dy_norm,
                               Dz = Dz_norm),
              image = V1_norm,
              actual_voxels = LL[[1]]$actual_voxels,
              voxel_lookup_table = LL[[1]]$voxel_lookup_table,
              candidates = my_candidates,
              patch_neighbours = LL[[1]]$patch_neighbours,
              similarities = my_similarities))

}
