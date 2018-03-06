build_dataset_for_ML <- function(problem_info,
                                 samples_per_subject = 1000,
                                 patch_size = 3,
                                 stride = 3,
                                 output_patch_size = 1,
                                 class_balance = TRUE,
                                 use_augmentation = TRUE) {

  if (!is.null(problem_info$train)) {

    # We have train and test sets

    train_set <- build_dataset(inputs = problem_info$train$x,
                               output = problem_info$train$y,
                               samples_per_subject = samples_per_subject,
                               patch_size = patch_size,
                               stride = stride,
                               output_patch_size = output_patch_size,
                               total_volumes = sum(problem_info$num_volumes),
                               class_balance = class_balance,
                               remap_classes = problem_info$remap_classes,
                               use_augmentation = use_augmentation)

    test_set <- build_dataset(inputs = problem_info$test$x,
                              output = problem_info$test$y,
                              samples_per_subject = samples_per_subject,
                              patch_size = patch_size,
                              stride = stride,
                              output_patch_size = output_patch_size,
                              total_volumes = sum(problem_info$num_volumes),
                              class_balance = class_balance,
                              remap_classes = problem_info$remap_classes,
                              use_augmentation = use_augmentation)

    return(list(train_set = train_set, test_set = test_set))

  } else {

    # Everything is training :-(

    dataset <- build_dataset(inputs = problem_info$inputs,
                             output = problem_info$outputs,
                             samples_per_subject = samples_per_subject,
                             patch_size = patch_size,
                             stride = stride,
                             output_patch_size = output_patch_size,
                             total_volumes = sum(problem_info$num_volumes),
                             class_balance = class_balance,
                             remap_classes = problem_info$remap_classes)

    return(list(dataset))

  }

}


build_dataset <- function(inputs,
                          output,
                          samples_per_subject = 1000,
                          patch_size = 3,
                          stride = 3,
                          output_patch_size = 1,
                          total_volumes = 1,
                          class_balance = TRUE,
                          remap_classes = NULL,
                          use_augmentation = TRUE) {

  subjects <- length(output)

  augmentation_times <- 1
  augment <- FALSE
  if (is.logical(use_augmentation) && use_augmentation) {

    augment <- TRUE

  }

  if (is.numeric(use_augmentation) && use_augmentation >= 1) {

    augmentation_times <- floor(use_augmentation)
    augment <- TRUE

  }

  # Dimensions of a set are:
  # rows: samples_per_subject * number_of subjects in the set * augmentation_times
  # columns: number of inputs of a subject (4d images as independent inputs) * patch_size ^ 3 + number of coordinates (3) +
  #           + output_patch_size ^ 3

  set <- matrix(0,
                nrow = samples_per_subject * subjects * augmentation_times,
                ncol = total_volumes * patch_size ^ 3 + 3 + output_patch_size ^ 3)

  if (!is.null(remap_classes))
    unique_labels <- unique(c(0, remap_classes$target, remap_classes$remaining))
  else {

    unique_labels <- 0
    class_balance <- FALSE

  }


  # Build train set
  for (s in seq(subjects)) {

    for (t in seq(augmentation_times)) {

      subject_data <- matrix(0, nrow = samples_per_subject * augmentation_times,
                             ncol = total_volumes * patch_size ^ 3 + 3 + output_patch_size ^ 3)

      if (augment) {

        M <- random_transformation_matrix()

      } else {

        M <- identity_matrix()

      }

      this_output <- read_nifti_to_array(output[s])

      if (!is.null(remap_classes)) {

        this_output <- map_ids_cpp(this_output, remap_classes)
        this_output <- apply_image_augmentation(this_output, M, type = "categorical")

      } else {

        this_output <- apply_image_augmentation(this_output, M, type = "continuous")

      }

      sample <- get_sample_indices(Vy = this_output,
                                   num_windows = samples_per_subject,
                                   batches_per_file = 1,
                                   stride = stride,
                                   mode = "sampling",
                                   class_balance = class_balance,
                                   unique_labels = unique_labels)

      sampling_indices <- sample$sampling_indices[which(sample$batch_idx == 1)]

      coords <- arrayInd(sampling_indices, .dim = dim(this_output))

      x <- coords[, 1] - 1
      y <- coords[, 2] - 1
      z <- coords[, 3] - 1

      my_outputs <- get_windows_at(V = this_output,
                                   width = output_patch_size,
                                   x = x, y = y, z = z)
      my_outputs <- my_outputs[, -c(1:3)]

      for (input in seq_along(inputs)) {

        # Read each input
        this_input <- read_nifti_to_array(inputs[[input]][s])
        this_input <- apply_image_augmentation(this_input, M, type = "continuous")

        X <- get_windows_at(V = this_input,
                            width = patch_size,
                            x = x, y = y, z = z)
        X <- X[, -c(1:3)]

        if (input == 1) {

          my_inputs <- X

        } else {

          my_inputs <- cbind(my_inputs, X)

        }

      }

      subject_data[seq((t - 1) * samples_per_subject + 1, t * samples_per_subject), ] <- cbind(coords - 1, my_inputs, my_outputs)

    }

    set[seq((s - 1) * samples_per_subject * augmentation_times + 1, s * samples_per_subject * augmentation_times), ] <- subject_data

  }

  return(set)

}
