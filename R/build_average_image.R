build_average_image <- function(image_list) {

  if (is.character(image_list)) {

    img <- read_nifti(image_list[1])
    n_images <- length(image_list)

    if (n_images > 1) {

      for (n in seq(2, n_images)) {

        img <- img + read_nifti(image_list[n])
      }

    }

    img <- img / n_images

    return(img)

  } else {

    n_images <- dim(image_list)[4]
    if (n_images > 1)
      img <- sum_4d(image_list)
    else
      img <- image_list

    return(img / n_images)

  }

}
