#' @title Read Multiple Files
#'
#' @description This function reads several NIfTI files.
#'
#' @param file_list    (list) Path of files to import
#'
#' @return A list with one entry for each file read.
#'
#' @export
#'
read_nifti_batch <- function(file_list) {

  return(lapply(file_list, read_nifti_to_array))

}
