sample_name_or_number_to_both <- function(sample, sample_names) {
  if (is.numeric(sample)) {
    sample_num <- sample
    if (any(sample_num < 0 || sample_num > length(sample_names))) {
      cli_abort("All samples should be between 1 and {length(sample_names)}")
    }
    sample_id <- sample_names[sample_num]
    logic <- rep(FALSE, length(sample_names))
    logic[sample_num] <- TRUE
  } else if (is.character(sample)) {
    sample_id <- sample
    sample_num <- match(sample, sample_names)
    missing_sample_names <- sample_id[which(is.na(sample_num))]
    if (length(missing_sample_names) > 0) {
      cli_abort("Missing sample names: {paste0(missing_sample_names, collapse = ', ')}")
    }
    logic <- rep(FALSE, length(sample_names))
    logic[sample_num] <- TRUE
  } else if (is.logical(sample)) {
    if (length(sample) != length(sample_names)) {
      cli_abort("When {.var sample} is a logical vector it should be of length {length(sample_names)} (the number of samples)")
    }
    logic <- sample
    logic[is.na(logic)] <- FALSE
    sample_id <- sample_names[logic]
    sample_num <- seq_along(sample_names)[logic]
  } else {
    cli_abort("{.var sample} should be a numeric vector, a character vector or a logical vector, not {class(sample)}")
  }

  list(
    idx = sample_num,
    name = sample_id,
    logical = logic
  )
}
