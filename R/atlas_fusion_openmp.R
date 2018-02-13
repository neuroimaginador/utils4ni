malf <- function(input_image,
                 template4D,
                 labels4D,
                 patch_size = 9,
                 search_size = 5,
                 k = 10,
                 stride = round((patch_size + 1) / 2),
                 max_iter = 10,
                 max_random_neighbours = 2,
                 label_ids = c(0, 1),
                 lambda = 0.7,
                 kernel_sigma = 1,
                 kernel_width = 3,
                 sigma2 = 0) {


  elapsed <- system.time(
    L <- obtain_candidates_similarities(image = input_image,
                                        templates = template4D,
                                        patch_size = patch_size,
                                        search_size = search_size,
                                        stride = stride,
                                        k = k,
                                        max_iter = max_iter,
                                        max_random_neighbours = max_random_neighbours)
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
                      k = k,
                      lambda = lambda,
                      sigma2 = sigma2,
                      match = L$similarities,
                      new_voting = new_voting,
                      ncores = 3)
  )
  cat("Elapsed in label fusion: ", elapsed[3], "\n")


  if (kernel_width > 0) {

    kernel_width <- 2 * floor(kernel_width / 2) + 1
    kernel <- gaussian_kernel(sigma = kernel_sigma, size = kernel_width)

    elapsed <- system.time(new_voting <- regularize(new_voting, kernel))
    cat("Elapsed in smoothing: ", elapsed[3], "\n")


  }

  seg <- defuzzify(new_voting) - 1

  return(seg)

}


obtain_candidates_similarities <- function(image,
                                           templates,
                                           patch_size = 9,
                                           search_size = 5,
                                           stride = round((patch_size + 1) / 2),
                                           k = 10,
                                           max_iter = 10,
                                           max_random_neighbours = 2) {


  voxel_lookup_table <- vector(mode = "integer", length = prod(dim(image)))

  actual_voxels <- count_elegible(image,
                                  patch_size = patch_size,
                                  search_size = search_size,
                                  stride = stride,
                                  voxel_lookup_table)

  kANN <- vector(mode = "integer", length = k * actual_voxels)
  similarities <- vector(mode = "numeric", length = k * actual_voxels)

  constrained_initialization_omp(input_image = image,
                                 template4D = templates,
                                 patch_size = patch_size,
                                 search_size = search_size,
                                 actual_voxels = actual_voxels,
                                 voxel_lookup_table = voxel_lookup_table,
                                 kANN = kANN,
                                 k = k,
                                 ncores = 3)



  patch_neighbours <- get_neighbours(array = image, width = patch_size)

  all_patches_similarity_omp(input_image = image,
                             template4D = templates,
                             k = k,
                             actual_voxels = actual_voxels,
                             voxel_lookup_table = voxel_lookup_table,
                             patch_neighbours = patch_neighbours,
                             kANN = kANN,
                             similarities = similarities,
                             ncores = 3)

  cat("Mean similarity = ", mean(similarities), "\n")

  direction <- -1
  for (i in seq(max_iter)) {

    propagation_step_omp(input_image = image,
                         template4D = templates,
                         actual_voxels = actual_voxels,
                         voxel_lookup_table = voxel_lookup_table,
                         patch_neighbours = patch_neighbours,
                         kANN = kANN,
                         k = k,
                         direction = direction,
                         patch_size = patch_size,
                         stride = stride,
                         similarities = similarities,
                         ncores = 3)

    cat("Mean similarity after PS = ", mean(similarities), "\n")

    direction <- -1 * direction

    constrained_random_search_omp(input_image = image,
                                  template4D = templates,
                                  actual_voxels = actual_voxels,
                                  voxel_lookup_table = voxel_lookup_table,
                                  kANN = kANN,
                                  k = k,
                                  patch_size = patch_size,
                                  patch_neighbours = patch_neighbours,
                                  search_size_max = search_size,
                                  similarities = similarities,
                                  max_random_neighbours = max_random_neighbours,
                                  ncores = 3)

    cat("Mean similarity after CRS = ", mean(similarities), "\n")

  }


  return(list(actual_voxels = actual_voxels,
              voxel_lookup_table = voxel_lookup_table,
              patch_neighbours = patch_neighbours,
              candidates = kANN,
              similarities = similarities))

}
