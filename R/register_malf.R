register_malf <- function(image, template_files, label_files, typeofTransform = "AffineFast", ...) {

  num_subjects <- length(template_files)

  system.time({

    templates <- array(0, dim = c(dim(image), num_subjects))
    labels <- array(0, dim = c(dim(image), num_subjects))

    if (!ANTsRCore::is.antsImage(image)) {

      img <- ANTsRCore::as.antsImage(image)

    } else {

     img <- image

    }

    for (i in seq(num_subjects)) {

      cat("Registering subject", i, "\n")

      tx <- ANTsRCore::antsRegistration(fixed = img,
                             moving = ANTsRCore::antsImageRead(template_files[i]),
                             verbose = FALSE,
                             typeofTransform = typeofTransform,
                             ...)

      templates[, , , i] <- tx$warpedmovout %>% as.array()

      labels[, , , i] <- ANTsRCore::antsApplyTransforms(fixed = img,
                                             moving = ANTsRCore::antsImageRead(label_files[i]),
                                             transformlist = tx$fwdtransforms,
                                             interpolator = "nearestNeighbor",
                                             verbose = FALSE) %>%
        as.array()

    }

  })

  return(list(templates = templates, labels = labels))

}
