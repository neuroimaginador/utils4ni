get_sample_indices <- function(Vy,
                               num_windows,
                               batches_per_file = 1,
                               stride = 5,
                               mode = c("sampling", "all"),
                               class_balance = TRUE,
                               unique_labels = c(0, 1),
                               verbose = FALSE) {

  mode <- mode[1]

  L <- expand.grid(seq(from = 1, to = dim(Vy)[1], by = stride),
                   seq(from = 1, to = dim(Vy)[2], by = stride),
                   seq(from = 1, to = dim(Vy)[3], by = stride))
  L <- as.matrix(L)

  all_idx <- L[, 1] + dim(Vy)[1] * (L[, 2] - 1) + dim(Vy)[1] * dim(Vy)[2] * (L[, 3] - 1)

  if (mode == "all") {

    sampling_indices <- all_idx

    if (verbose)
      message("Number of actual windows: ", length(sampling_indices)) # nocov

    num_batches <- ceiling(length(sampling_indices) / num_windows)
    max_epochs <- min(c(num_batches, batches_per_file))

    if (max_epochs > 1) {

      batch_idx <- rep(seq(max_epochs), each = num_windows)
      batch_idx <- batch_idx[seq_along(sampling_indices)]

    } else
      batch_idx <- rep(1, times = length(sampling_indices))

  } else {

    sampling_indices <- sample(all_idx, length(all_idx))

    if (class_balance & (length(unique_labels) > 1)) {

      # Vy <- map_ids(image = Vy, remap_classes = config$remap_classes)
      # unique_labels <- unique(c(0, config$remap_classes$target, config$remap_classes$remaining))

      balanced_classes <- sample(unique_labels, size = length(all_idx), replace = TRUE)
      sampling_indices <- rep(0, length(balanced_classes))

      for (class in unique_labels) {

        idx_for_class <- which(balanced_classes == class)
        idx_in_img <- which(Vy == class)

        if (length(idx_in_img) > 0) {
          idx <- sample(idx_in_img,
                        size = length(idx_for_class),
                        replace = length(idx_for_class) > length(idx_in_img))
          sampling_indices[idx_for_class] <- idx

        }

      }

      sampling_indices <- sampling_indices[sampling_indices > 0]

      if (length(sampling_indices) < length(all_idx))
        sampling_indices <- sample(sampling_indices,
                                   size = length(all_idx),
                                   replace = TRUE)


    }

    num_batches <- ceiling(length(sampling_indices) / num_windows)
    max_epochs <- min(c(num_batches, batches_per_file))

    batch_idx <- rep(seq(max_epochs), each = num_windows)

  }

  return(list(sampling_indices = sampling_indices,
              num_batches = num_batches,
              max_epochs = max_epochs,
              batch_idx = batch_idx))

}
