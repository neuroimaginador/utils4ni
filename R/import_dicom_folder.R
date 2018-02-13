import_dicom_folder <- function(folder) {

  require(divest)

  image_list <- readDicom(path = folder, interactive = FALSE)

  if (length(image_list) == 0) {

    cmd <- paste0("find ", folder, ' -name "*" -exec gdcmconv -w {} {} \\;')
    system(cmd)

    image_list <- readDicom(path = folder, interactive = FALSE)

  }

  if (length(image_list) > 0)
    return(image_list)
  else
    stop("No available DICOM in such folder.")

}
