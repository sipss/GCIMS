#' Sample names
#'
#' @param object GCIMSDataset object
#' @param value A character vector of length the number of samples with the sample names
#' @return The [GCIMSDataset] object
#' @importMethodsFrom Biobase sampleNames
#' @importMethodsFrom Biobase "sampleNames<-"
setMethod("sampleNames", "GCIMSDataset", function(object) object$pData$SampleID)

#' @describeIn sampleNames-GCIMSDataset-method Sample names
setReplaceMethod("sampleNames", "GCIMSDataset", function(object, value) {
  if (nrow(object$pData) != length(value)) {
    cli_abort(
      c(
        "Invalid sample names",
        "x" = "The number of sample names given ({length(value)}) != Number of samples ({nrow(object$pData)})"
      )
    )
  }
  if (anyNA(value) || anyDuplicated(value)) {
    cli_abort("Sample names must be unique and not missing")
  }
  rename_intermediate_files(object, value)
  object$pData$SampleID <- value
  object
})

rename_intermediate_files <- function(object, new_names) {
  current_dir <- object$currentHashedDir()
  if (is.null(current_dir)) {
    return()
  }
  current_dir_old <- paste0(current_dir, "_old")
  if (!dir.exists(current_dir)) {
    return()
  }
  file.rename(current_dir, current_dir_old)
  dir.create(current_dir, showWarnings = FALSE, recursive = TRUE)
  old_names <- names(object)
  for (i in seq_along(new_names)) {
    old_name <- paste0(old_names[i], ".rds")
    new_name <- paste0(new_names[i], ".rds")
    if (file.exists(file.path(current_dir_old, old_name))) {
      file.rename(
        file.path(current_dir_old, old_name),
        file.path(current_dir, new_name)
      )
    }
  }
  unlink(current_dir_old, recursive = TRUE)
  invisible(NULL)
}
