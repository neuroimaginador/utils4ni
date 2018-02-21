obtain_candidates_files <- function(image,
                                    mask = NULL,
                                    templates,
                                    patch_size = 9,
                                    search_size = 5,
                                    stride = round((patch_size + 1) / 2),
                                    max_iter = 10,
                                    max_random_neighbours = 2,
                                    ncores = 2) {

  n_templates <- length(templates)

  voxel_lookup_table <- vector(mode = "integer", length = prod(dim(image)))

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

  k <- n_templates
  kANN <- vector(mode = "integer", length = k * actual_voxels)
  similarities <- vector(mode = "numeric", length = k * actual_voxels)

  patch_neighbours <- get_neighbours(array = image, width = patch_size)

  for (k1 in seq_along(templates)) {

    template <- read_nifti(templates[k1])
    dim(template) <- c(dim(template), 1)
    kANN_temp <- kANN[(k1 - 1) * actual_voxels + seq(actual_voxels)]
    simil_temp <- similarities[(k1 - 1) * actual_voxels + seq(actual_voxels)]

    constrained_initialization_omp(input_image = image,
                                   template4D = template,
                                   patch_size = patch_size,
                                   search_size = search_size,
                                   actual_voxels = actual_voxels,
                                   voxel_lookup_table = voxel_lookup_table,
                                   kANN = kANN_temp,
                                   ncores = ncores)

    all_patches_similarity_omp(input_image = image,
                               template4D = template,
                               actual_voxels = actual_voxels,
                               voxel_lookup_table = voxel_lookup_table,
                               patch_neighbours = patch_neighbours,
                               kANN = kANN_temp,
                               similarities = simil_temp,
                               ncores = ncores)

    cat("Mean similarity = ", mean(simil_temp), "\n")

    direction <- -1

    for (i in seq(max_iter)) {

      propagation_step_omp(input_image = image,
                           template4D = template,
                           actual_voxels = actual_voxels,
                           voxel_lookup_table = voxel_lookup_table,
                           patch_neighbours = patch_neighbours,
                           kANN = kANN_temp,
                           direction = direction,
                           patch_size = patch_size,
                           stride = stride,
                           similarities = simil_temp,
                           ncores = ncores)

      cat("Mean similarity after PS- = ", mean(simil_temp), "\n")

      direction <- -1 * direction

      constrained_random_search_omp(input_image = image,
                                    template4D = template,
                                    actual_voxels = actual_voxels,
                                    voxel_lookup_table = voxel_lookup_table,
                                    kANN = kANN_temp,
                                    patch_size = patch_size,
                                    patch_neighbours = patch_neighbours,
                                    search_size_max = search_size,
                                    similarities = simil_temp,
                                    max_random_neighbours = max_random_neighbours,
                                    ncores = ncores)

      cat("Mean similarity after CRS = ", mean(simil_temp), "\n")

    }

    kANN[(k1 - 1) * actual_voxels + seq(actual_voxels)] <- kANN_temp
    similarities[(k1 - 1) * actual_voxels + seq(actual_voxels)] <- simil_temp

  }



  return(list(actual_voxels = actual_voxels,
              voxel_lookup_table = voxel_lookup_table,
              patch_neighbours = patch_neighbours,
              candidates = kANN,
              similarities = similarities))

}
