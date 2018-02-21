malf_memory <- function(input_image,
                        templates,
                        labels,
                        patch_size = 9,
                        search_size = 5,
                        stride = round((patch_size + 1) / 2),
                        k = 10,
                        max_iter = 10,
                        max_random_neighbours = 2,
                        label_ids = c(0, 1),
                        lambda = 0.7,
                        kernel_sigma = 1,
                        kernel_width = 3,
                        sigma2 = 0,
                        ncores = parallel::detectCores() - 1) {

  cat("Running MALF in", ncores, "cores,\n")

  elapsed <- system.time(
    L <- obtain_candidates(input_image = input_image,
                           templates = templates,
                           patch_size = patch_size,
                           search_size = search_size,
                           stride = stride,
                           k = k,
                           max_iter = max_iter,
                           max_random_neighbours = max_random_neighbours,
                           ncores = ncores)
  )

  cat("Elapsed in obtaining candidates: ", elapsed[3], "\n")

  label_ids <- label_ids %>% as.integer()
  new_voting <- array(0, dim = c(dim(input_image), length(label_ids)))

  elapsed <- system.time(
    label_fusion2_omp(labels4D = labels,
                      actual_voxels = L$actual_voxels %>% as.integer(),
                      voxel_lookup_table = L$voxel_lookup_table %>% as.integer(),
                      k = k,
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

  seg <- defuzzify(new_voting) - 1

  return(seg)

}


obtain_candidates_similarities_memory <- function(image,
                                                  template_files,
                                                  patch_size = 9,
                                                  search_size = 5,
                                                  stride = round((patch_size + 1) / 2),
                                                  max_iter = 10,
                                                  max_random_neighbours = 2,
                                                  ncores = 2) {

  # image <- map_images(source = image, target = templates[, , , 1])

  voxel_lookup_table <- vector(mode = "integer", length = prod(dim(image)))

  cat("P = ", prod(dim(image)), "\n")

  cat("Actual voxels\n")
  actual_voxels <- count_elegible(image,
                                  patch_size = patch_size,
                                  search_size = search_size,
                                  stride = stride,
                                  voxel_lookup_table = voxel_lookup_table)
  print(actual_voxels)

  k <- length(template_files)

  kANN <- vector(mode = "integer", length = k * actual_voxels)
  similarities <- vector(mode = "numeric", length = k * actual_voxels)

  print(k * actual_voxels)
  print(length(kANN))

  Sys.sleep(5)

  # stop("ERROR GORDO")

  cat("Constrained init\n")
  constrained_init_memory(input_image = image,
                          n_templates = length(template_files),
                          patch_size = patch_size,
                          search_size = search_size,
                          actual_voxels = actual_voxels,
                          voxel_lookup_table = voxel_lookup_table,
                          kANN = kANN,
                          ncores = ncores)

  cat("Neighbours\n")

  patch_neighbours <- get_neighbours(array = image, width = patch_size)

  cat("Patches similarity\n")
  patches_similarity_memory(input_image = image,
                            template_filenames = template_files,
                            actual_voxels = actual_voxels,
                            voxel_lookup_table = voxel_lookup_table,
                            patch_neighbours = patch_neighbours,
                            kANN = kANN,
                            similarities = similarities,
                            ncores = ncores)

  cat("Mean similarity = ", mean(similarities), "\n")

  direction <- -1
  for (i in seq(max_iter)) {

    propagation_step_memory(input_image = image,
                            template_filenames = template_files,
                            actual_voxels = actual_voxels,
                            voxel_lookup_table = voxel_lookup_table,
                            patch_neighbours = patch_neighbours,
                            kANN = kANN,
                            direction = direction,
                            patch_size = patch_size,
                            stride = stride,
                            similarities = similarities,
                            ncores = ncores)

    cat("Mean similarity after PS = ", mean(similarities), "\n")

    direction <- -1 * direction

    constrained_random_search_memory(input_image = image,
                                     template_filenames = template_files,
                                     actual_voxels = actual_voxels,
                                     voxel_lookup_table = voxel_lookup_table,
                                     kANN = kANN,
                                     patch_size = patch_size,
                                     patch_neighbours = patch_neighbours,
                                     search_size_max = search_size,
                                     similarities = similarities,
                                     max_random_neighbours = max_random_neighbours,
                                     ncores = ncores)

    cat("Mean similarity after CRS = ", mean(similarities), "\n")

  }


  return(list(actual_voxels = actual_voxels,
              voxel_lookup_table = voxel_lookup_table,
              patch_neighbours = patch_neighbours,
              candidates = kANN,
              similarities = similarities))

}
