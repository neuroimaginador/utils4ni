#' Get (Linear) Coordinates of Neighborhood Voxels
#'
#' @param x       (Array) Reference image on where to compute neighborhood
#' @param width   (Integer) Width of the neighborhood
#' @param type    (Character string) A string giving the type of shape to produce. In one dimension, these shapes are all equivalent.
#'
#' @return A list with two components: \code{neighbours}, the list of indices of the neighbours, and \code{distances}, The distances from the center voxel to each of the neighbours.
#'
get_neighborhood <- function(x, width, type = c("box", "disc", "diamond")) {
  
  # Dependencies
  require(mmand)
  require(oro.nifti)
  require(ANTsR)
  
  # Initialization of variables
  dims <- dim(x)
  voxel_dims <- rep(1, length(dims))
  if (is.nifti(x)) voxel_dims <- pixdim(x)[2:4]
  if (is.antsImage(x)) voxel_dims <- antsGetSpacing(x)
  
  # Kernel to use and neighbours
  kernel <- shapeKernel(width = rep(width, length.out = length(dims)), 
                        type = type)
  
  neighbours <- neighbourhood(x, width = width)$offset[as.vector(kernel) > 0]
  
  # Linear indices of the neighbours
  res <- arrayInd(which(kernel > 0), .dim = dim(kernel))
  
  # Computation of distances
  idx_center <- which(neighbours == 0)
  
  for (i in 1:length(dims)) {
    
    res[, i] <- (res[, i] - res[idx_center, i]) * voxel_dims[i]
    
  }
  
  distances <- apply(res, 1, FUN = function(row) sum(row ^ 2))
  
  return(list(neighbours = neighbours, distances = distances))
  
}
