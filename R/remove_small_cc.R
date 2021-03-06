remove_small_cc <- function(labelling, pctg = 0.25) {


  cc <- connected_components(labelling)
  labels <- sort(unique(as.vector(cc)))
  labels <- labels[labels > 0]

  remap_classes <- list(source = labels, target = seq_along(labels))

  cc <- map_ids_cpp(cc, remap_classes = remap_classes, invert = FALSE)

  # confusion_matrix <- table(labelling, cc)
  cm <- confusion_matrix(labelling, cc)
  cm <- cm[-1, -1]

  if (is.vector(cm)) {

    idx_below <- which_below(cm, pctg)

  } else {

    idx_below <- unlist(lapply(1:nrow(cm),
                               function(row) {

                                 v <- cm[row, ]
                                 which_below(v, pctg)

                               }))

  }

  remap_classes <- list(source = idx_below,
                        target = rep(1, length.out = length(idx_below)),
                        remaining = 0)

  out_mask <- map_ids_cpp(cc, remap_classes = remap_classes)

  labelling[out_mask > 0] <- 0

  return(labelling)

}

which_below <- function(v, pctg) {

  max_component_size <- max(v)
  threshold <- pctg * max_component_size
  which(v > 0 & v <= threshold)

}
