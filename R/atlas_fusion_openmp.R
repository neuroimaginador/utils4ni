malf <- function(input_image,
                 mask = NULL,
                 template4D,
                 labels4D,
                 patch_size = 9,
                 search_size = 5,
                 stride = round((patch_size + 1) / 2),
                 max_iter = 10,
                 max_random_neighbours = 2,
                 label_ids = c(0, 1),
                 lambda = 0.7,
                 kernel_sigma = 1,
                 kernel_width = 3,
                 sigma2 = 0,
                 method = c("mi", "ssd"),
                 early_stopping = list(),
                 return_memberships = FALSE,
                 ncores = parallel::detectCores() - 1) {

  cat("Running MALF in", ncores, "cores\n")
  cat("Using", method[1], "as similarity.\n")

  elapsed <- system.time(
    L <- obtain_candidates_similarities(image = input_image,
                                        mask = mask,
                                        templates = template4D,
                                        patch_size = patch_size,
                                        search_size = search_size,
                                        stride = stride,
                                        max_iter = max_iter,
                                        method = method,
                                        max_random_neighbours = max_random_neighbours,
                                        early_stopping = early_stopping,
                                        ncores = ncores)
  )

  cat("Elapsed in obtaining candidates: ", elapsed[3], "\n")

  label_ids <- label_ids %>% as.integer()

  new_voting <- array(0, dim = c(dim(input_image), length(label_ids)))

  elapsed <- system.time(
    label_fusion2_omp(labels4D = labels4D,
                      actual_voxels = L$actual_voxels %>% as.integer(),
                      voxel_lookup_table = L$voxel_lookup_table %>% as.integer(),
                      label_ids = label_ids,
                      kANN = L$candidates %>% as.integer(),
                      patch_neighbours = L$patch_neighbours %>% as.integer(),
                      lambda = lambda,
                      sigma2 = sigma2,
                      match = L$similarities,
                      new_voting = new_voting,
                      ncores = ncores)
  )
  cat("Elapsed in label fusion: ", elapsed[3], "\n")


  if (kernel_width > 0) {

    kernel_width <- 2 * floor(kernel_width / 2) + 1
    kernel <- gaussian_kernel(sigma = kernel_sigma, size = kernel_width)

    elapsed <- system.time(new_voting <- regularize(new_voting, kernel, ncores = ncores))
    cat("Elapsed in smoothing: ", elapsed[3], "\n")


  }

  if (return_memberships) {

    return(new_voting)

  } else {

    seg <- defuzzify(new_voting) - 1

    return(seg)

  }

}

fast_malf <- function(input_image,
                      mask = NULL,
                      template4D,
                      labels4D,
                      patch_size = 9,
                      search_size = 5,
                      stride = round((patch_size + 1) / 2),
                      max_iter = 10,
                      max_random_neighbours = 2,
                      label_ids = c(0, 1),
                      lambda = 0.7,
                      kernel_sigma = 1,
                      kernel_width = 3,
                      sigma2 = 0,
                      method = c("mi", "ssd"),
                      ncores = parallel::detectCores() - 1) {

  cat("Running MALF in", ncores, "cores\n")
  cat("Using", method[1], "as similarity.\n")

  elapsed <- system.time(
    L <- obtain_candidates_similarities(image = input_image,
                                        mask = mask,
                                        templates = template4D,
                                        patch_size = patch_size,
                                        search_size = search_size,
                                        stride = stride,
                                        max_iter = max_iter,
                                        method = method,
                                        max_random_neighbours = max_random_neighbours,
                                        ncores = ncores)
  )

  cat("Elapsed in obtaining candidates: ", elapsed[3], "\n")

  label_ids <- label_ids %>% as.integer()

  cat("Fast Mode\n")

  new_sim <- array(0, dim = dim(labels4D))
  dim(new_sim) <- dim(labels4D)

  new_voting <- new_sim %>% as.integer()
  dim(new_voting) <- dim(labels4D)

  elapsed <- system.time( {
    label_fusion_omp_fast(labels4D = labels4D,
                          actual_voxels = L$actual_voxels %>% as.integer(),
                          voxel_lookup_table = L$voxel_lookup_table %>% as.integer(),
                          label_ids = label_ids,
                          kANN = L$candidates %>% as.integer(),
                          patch_neighbours = L$patch_neighbours %>% as.integer(),
                          lambda = lambda,
                          sigma2 = sigma2,
                          match = L$similarities,
                          new_voting = new_voting,
                          new_sim = new_sim,
                          ncores = ncores)

    cat("Label Fusion\n")

    gc()

    if (kernel_width > 0) {

      cat("Regularization\n")

      kernel_width <- 2 * floor(kernel_width / 2) + 1

      kernel <- gaussian_kernel(sigma = kernel_sigma,
                                dim = 3,
                                size = kernel_width,
                                normalised = TRUE)

      new_sim <- regularize(image = new_sim,
                            kernel = kernel,
                            ncores = ncores)

    }

    result <- array(0, dim = dim(input_image)) %>% as.integer()

    label_fusion_mode(my_labels = new_voting,
                      my_similarities = new_sim,
                      result = result,
                      ncores = ncores)

    dim(result) <- dim(input_image)

    cat("Mode\n")

  })

  cat("Elapsed in label fusion: ", elapsed[3], "\n")

  return(result)

}



obtain_candidates_similarities <- function(image,
                                           mask = NULL,
                                           templates,
                                           patch_size = 9,
                                           search_size = 5,
                                           stride = round((patch_size + 1) / 2),
                                           max_iter = 10,
                                           max_random_neighbours = 2,
                                           crs = TRUE,
                                           early_stopping = list(),
                                           scheme = 0,
                                           method = c("mi", "ssd"),
                                           ncores = 2) {

  # image <- map_images(source = image, target = templates[, , , 1])
  #
  method <- method[1]
  my_method <- switch(method,
                      "mi" = 0,
                      "ssd" = 1)

  voxel_lookup_table <- vector(mode = "integer", length = prod(dim(image)))
  cat("Init\n")

  if (is.null(mask)) {

    actual_voxels <- count_elegible(image,
                                    patch_size = patch_size,
                                    search_size = search_size,
                                    stride = stride,
                                    voxel_lookup_table)

  } else {

    actual_voxels <- count_elegible_masked(image,
                                           mask,
                                           patch_size = patch_size,
                                           search_size = search_size,
                                           stride = stride,
                                           voxel_lookup_table)

  }

  k <- dim(templates)[4]
  kANN <- vector(mode = "integer", length = k * actual_voxels)
  similarities <- vector(mode = "numeric", length = k * actual_voxels)

  cat("k = ", k, "\n")
  cat("actual_voxels = ", actual_voxels, "\n")

  cat("Constrained init\n")

  patch_neighbours <- get_neighbours(array = image, width = patch_size)

  constrained_initialization_omp(input_image = image,
                                 template4D = templates,
                                 patch_size = patch_size,
                                 search_size = search_size,
                                 actual_voxels = actual_voxels,
                                 voxel_lookup_table = voxel_lookup_table,
                                 kANN = kANN,
                                 ncores = ncores)

  cat("Iterations\n")

  for (sc in scheme) {

    cat("Initial scheme : ", sc, "\n")

    if (sc > 0) {

      kernel <- gaussian_kernel(sigma = sc, size = 5)
      this_image <- regularize(image, kernel, ncores = ncores)
      these_templates <- regularize(templates, kernel)

    } else {

      this_image <- image
      these_templates <- templates

    }

    all_patches_similarity_omp(input_image = this_image,
                               template4D = these_templates,
                               actual_voxels = actual_voxels,
                               voxel_lookup_table = voxel_lookup_table,
                               patch_neighbours = patch_neighbours,
                               kANN = kANN,
                               similarities = similarities,
                               method = my_method,
                               ncores = ncores)

    previous_similarity <- mean(similarities)
    cat("Mean similarity = ", previous_similarity, "\n")

    direction <- -1

    for (i in seq(max_iter)) {

      propagation_step_omp(input_image = this_image,
                           template4D = these_templates,
                           actual_voxels = actual_voxels,
                           voxel_lookup_table = voxel_lookup_table,
                           patch_neighbours = patch_neighbours,
                           kANN = kANN,
                           direction = direction,
                           patch_size = patch_size,
                           stride = stride,
                           similarities = similarities,
                           method = my_method,
                           ncores = ncores)

      cat("Mean similarity after PS- = ", mean(similarities), "\n")

      direction <- -1 * direction

      # if (!crs) {
      propagation_step_omp(input_image = this_image,
                           template4D = these_templates,
                           actual_voxels = actual_voxels,
                           voxel_lookup_table = voxel_lookup_table,
                           patch_neighbours = patch_neighbours,
                           kANN = kANN,
                           direction = direction,
                           patch_size = patch_size,
                           stride = stride,
                           similarities = similarities,
                           method = my_method,
                           ncores = ncores)


      cat("Mean similarity after PS+ = ", mean(similarities), "\n")
      direction <- -1 * direction

      # }

      if (crs) {

        constrained_random_search_omp(input_image = this_image,
                                      template4D = these_templates,
                                      actual_voxels = actual_voxels,
                                      voxel_lookup_table = voxel_lookup_table,
                                      kANN = kANN,
                                      patch_size = patch_size,
                                      patch_neighbours = patch_neighbours,
                                      search_size_max = search_size,
                                      similarities = similarities,
                                      max_random_neighbours = max_random_neighbours,
                                      method = my_method,
                                      ncores = ncores)

        cat("Mean similarity after CRS = ", mean(similarities), "\n")

      }

      current_similarity <- mean(similarities)
      difference <- abs((previous_similarity - current_similarity) / previous_similarity)
      # cat("Difference = ", difference, "\n")
      # str(early_stopping)
      if ("tol" %in% names(early_stopping)) {

        tol <- as.numeric(early_stopping$tol)

        if (difference < tol) {

          cat("Reached maximum allowed tolerance.\n")
          break

        }

      }
      previous_similarity <- current_similarity

    }


  }


  return(list(actual_voxels = actual_voxels,
              voxel_lookup_table = voxel_lookup_table,
              patch_neighbours = patch_neighbours,
              candidates = kANN,
              similarities = similarities))

}

# constrained_init <- function(patch_neighbours, voxel_lookup_table, actual_voxels, n_templates) {
#
#   v <- which(voxel_lookup_table >= 0) - 1
#   idx <- voxel_lookup_table[v + 1]
#
#   displ <- sample(patch_neighbours, replace = TRUE, size = n_templates * actual_voxels)
#
#   res <- rep(v, times = n_templates) + displ
#   res[res < 0] <- 0
#   res[res >= length(voxel_lookup_table)] <- length(voxel_lookup_table) - 1
#
#   return(res)
#
# }
