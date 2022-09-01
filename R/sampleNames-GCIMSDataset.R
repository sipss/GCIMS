#' Sample names
#'
#' @param object GCIMSDataset object
#' @param value A character vector of length the number of samples with the sample names
#' @return The [GCIMSDataset] object
#' @importMethodsFrom Biobase sampleNames
#' @importMethodsFrom Biobase "sampleNames<-"
setMethod("sampleNames", "GCIMSDataset", function(object) object@envir$pData$SampleID)

#' @describeIn sampleNames-GCIMSDataset-method Sample names
setReplaceMethod("sampleNames", "GCIMSDataset", function(object, value) {
  if (nrow(object@envir$pData) != length(value)) {
    rlang::abort(
      message = c(
        "Invalid sample names",
        glue("The length of the given sample names ({length(value)}) should be equal to the number of samples ({nrow(object@envir$pData)})")
      )
    )
  }
  if (anyDuplicated(value)) {
    rlang::abort("Sample names must be unique")
  }
  rename_intermediate_files(object, value)
  object@envir$pData$SampleID <- value
  object
})

rename_intermediate_files <- function(object, new_names) {
  current_dir <- CurrentHashedDir(object)
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
