move_samples_current_to_prev <- function(object) {
  prev_dir <- file.path(object@envir$scratch_dir, "samples_prev")
  current_dir <- file.path(object@envir$scratch_dir, "samples_now")


  if (dir.exists(prev_dir)) {
    unlink(prev_dir, recursive = TRUE)
  }
  if (dir.exists(current_dir)) {
    file.rename(current_dir, prev_dir)
  }
  dir.create(current_dir, showWarnings = FALSE, recursive = TRUE)
  invisible(NULL)
}


realize_one_sample <- function(sample_name, prev_dir, current_dir, delayed_ops, base_dir, orig_filenames) {
  out <- vector("list", length = length(delayed_ops))
  needs_re_saving <- FALSE

  f <- paste0(sample_name, ".rds")
  sample_fn_prev <- file.path(prev_dir, f)
  sample_fn_curr <- file.path(current_dir, f)
  if (!file.exists(sample_fn_prev)) {
    # If the file does not exist and we are not reading the samples, where are we?
    if (delayed_ops[[1]]@name != "read_sample") {
      rlang::abort(
        message = c(
          "UnexpectedError",
          "x" = glue("Sample {sample_name} should have a file at {sample_fn_prev} but it was not found"),
          "x" = glue("Or the first action on the dataset should have been read_samples but it was not"),
          "i" = "Either case, this is a problem. You can try deleting the scratch directory and restart again"
        )
      )
    }
    # Read the input file:
    gcimssample <- file.path(base_dir, orig_filenames[sample_name])
  } else {
    # Read the file from the previous output:
    if (delayed_ops[[1]]@name == "read_sample") {
      rlang::abort(
        message = c(
          "x" = "The scratch dir already has samples loaded, but we are reading the samples",
          "i" = "Please delete the scratch dir and try again"
        )
      )
    }
    gcimssample <- readRDS(sample_fn_prev)
    if (!identical(gcimssample@description, sample_name)) {
      # sampleNames were probably updated, files renamed.
      gcimssample@description <- sample_name
      needs_re_saving <- TRUE
    }
  }
  for (i in seq_along(delayed_ops)) {
    result <- apply_op_to_sample(delayed_ops[[i]], gcimssample)
    gcimssample <- result$sample
    if (!identical(gcimssample@description, sample_name)) {
      # sampleNames were probably updated, files renamed.
      gcimssample@description <- sample_name
      needs_re_saving <- TRUE
    }
    if (!is.null(result$extracted_obj)) {
      out[[i]] <- result$extracted_obj
    }
    if (result$needs_resaving) {
      needs_re_saving <- TRUE
    }
  }
  # saveRDS, or just copy
  if (needs_re_saving) {
    saveRDS(gcimssample, sample_fn_curr)
  } else {
    if (!file.copy(
      sample_fn_prev,
      sample_fn_curr
    )) {
      rlang::abort(glue("Could not copy {sample_fn_prev} to {sample_fn_curr}."))
    }
  }
  out
}


#' Runs all delayed operations on the object
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @return The same [GCIMSDataset] object, without pending operations
#' @export
#' @examples
#' base_dir <- system.file("extdata", "sample_formats", package = "GCIMS")
#' annot <- data.frame(SampleID = "Sample1", FileName = "small.mea.gz")
#' dataset <- GCIMSDataset(annot, base_dir)
#' print(dataset)
#' realize(dataset)
#' print(dataset)
#'
realize <- function(object) {
  if (!hasDelayedOps(object)) {
    return(object)
  }
  move_samples_current_to_prev(object)
  prev_dir <- file.path(object@envir$scratch_dir, "samples_prev")
  current_dir <- file.path(object@envir$scratch_dir, "samples_now")

  sample_names <- sampleNames(object)
  names(sample_names) <- sample_names

  pdata <- Biobase::pData(object)
  orig_filenames <- pdata$FileName
  names(orig_filenames) <- sample_names

  extracted_results <- BiocParallel::bplapply(
    X = sample_names,
    FUN = realize_one_sample,
    prev_dir = prev_dir,
    current_dir = current_dir,
    delayed_ops = object@envir$delayed_ops,
    base_dir = object@envir$base_dir,
    orig_filenames = orig_filenames
  )

  # Apply to the dataset object
  for (i in seq_along(object@envir$delayed_ops)) {
    delayed_op <- object@envir$delayed_ops[[i]]
    # Extract i-th result for all samples:
    extracted_result <- purrr::map(extracted_results, i)
    # Let the delayed operation aggregate the extracted results and save them in object
    object <- aggregate_result(delayed_op, extracted_result, object)
  }
  # Make delayed_ops history and remove them from pending
  object@envir$previous_ops <- c(object@envir$previous_ops, object@envir$delayed_ops)
  object@envir$delayed_ops <- list()
  invisible(object)
}
