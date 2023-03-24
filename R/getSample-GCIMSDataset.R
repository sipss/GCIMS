#' Get a sample from a GCIMSDataset object
#'
#' @param sample Either an integer (sample index) or a string (sample name)
#' @return The corresponding [GCIMSSample]
setMethod(
  "getSample",
  "GCIMSDataset",
  function(object, sample) {
    object$getSample(sample)
  }
)


sample_name_or_number_to_both <- function(sample, sample_names) {
  if (is.numeric(sample)) {
    sample_num <- sample
    if (any(sample_num < 0 || sample_num > length(sample_names))) {
      cli_abort("All samples should be between 1 and {length(sample_names)}")
    }
    sample_id <- sample_names[sample_num]
  } else if (is.character(sample)) {
    sample_id <- sample
    sample_num <- match(sample, sample_names)
    missing_sample_names <- sample_id[which(is.na(sample_num))]
    if (length(missing_sample_names) > 0) {
      cli_abort("Missing sample names: {paste0(missing_sample_names, collapse = ', ')}")
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
