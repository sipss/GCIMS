#' @describeIn GCIMSDataset class
#' @param sample Either an integer (sample index) or a string (sample name)
#' @return The corresponding [GCIMSSample]
setMethod(
  "getGCIMSSample",
  "GCIMSDataset",
  function(object, sample) {
    object <- realize(object)
    if (is.numeric(sample)) {
      if (sample > length(sampleNames(object))) {
        stop("sample exceeds the number of samples")
      }
      sampleid <- sampleNames(object)[sample]
    } else if (is.character(sample)) {
      if (!sample %in% sampleNames(object)) {
        stop("sample not in sampleNames")
      }
      sampleid <- sample
    }
    filename <- paste0(sampleid, ".rds")
    sample_file <- file.path(object@envir$scratch_dir, "samples_now", filename)
    if (!file.exists(sample_file)) {
      rlang::abort(glue("File not found: {sample_file} should have been created"))
    }
    gcimssample <- readRDS(sample_file)
    if (!methods::is(gcimssample, "GCIMSSample")) {
      rlang::abort("Expected a GCIMSSample object, but it was not found")
    }
    gcimssample
  }
)

