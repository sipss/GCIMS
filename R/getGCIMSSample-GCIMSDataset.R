#' @describeIn GCIMSDataset class
#' @param sample Either an integer (sample index) or a string (sample name)
#' @return The corresponding [GCIMSSample]
setMethod(
  "getGCIMSSample",
  "GCIMSDataset",
  function(object, sample) {
    object <- realize(object)
    sample_id_num <- sample_name_or_number_to_both(sample, sampleNames(object))
    if (object@envir$on_ram) {
      gcimssample <- object@envir$samples[[sample_id_num$idx]]
    } else {
      filename <- paste0(sample_id_num$name, ".rds")
      current_intermediate_dir <- CurrentHashedDir(object)
      sample_file <- file.path(current_intermediate_dir, filename)
      if (!file.exists(sample_file)) {
        abort(glue("File not found: {sample_file} should have been created"))
      }
      gcimssample <- readRDS(sample_file)
    }
    if (!methods::is(gcimssample, "GCIMSSample")) {
      abort("Expected a GCIMSSample object, but it was not found")
    }
    gcimssample
  }
)


sample_name_or_number_to_both <- function(sample, sample_names) {
  if (is.numeric(sample)) {
    sample_num <- sample
    if (any(sample_num < 0 || sample_num > length(sample_names))) {
      abort(glue("All samples should be between 1 and {length(sample_names)}"))
    }
    sample_id <- sample_names[sample_num]
  } else if (is.character(sample)) {
    sample_id <- sample
    sample_num <- match(sample, sample_names)
    missing_sample_names <- sample_id[which(is.na(sample_num))]
    if (length(missing_sample_names) > 0) {
      abort(glue("Missing sample names: {paste0(missing_sample_names, collapse = ', ')}"))
    }
  }
  logic <- rep(FALSE, length(sample_names))
  logic[sample_num] <- TRUE

  list(
    idx = sample_num,
    name = sample_id,
    logical = logic
  )
}
