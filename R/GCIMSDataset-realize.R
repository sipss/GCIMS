
realize_one_sample <- function(sample_name, curr_dir, next_dir, delayed_ops, base_dir, orig_filenames) {
  if (is.null(next_dir)) {
    rlang::abort(message = c("UnexpectedError", "x" = "next_dir should not have been NULL"))
  }

  out <- vector("list", length = length(delayed_ops))
  needs_re_saving <- FALSE

  f <- paste0(sample_name, ".rds")
  sample_fn_next <- file.path(next_dir, f)

  if (is.null(curr_dir) && delayed_ops[[1]]@name != "read_sample") {
    rlang::abort(
      message = c(
        "UnexpectedError",
        "x" = glue("The first operation should have been named read_sample instead of {delayed_ops[[1]]@name}"),
        "i" = "This is an unexpected problem. You can try deleting the scratch directory and restart again"
      )
    )
  }

  if (is.null(curr_dir)) {
    # we are going to read_sample
    gcimssample <- file.path(base_dir, orig_filenames[sample_name])
    sample_fn_curr <- NULL
  } else {
    # We load the file
    sample_fn_curr <- file.path(curr_dir, f)
    if (!file.exists(sample_fn_curr)) {
      rlang::abort(
        message = c(
          "UnexpectedError",
          "x" = glue("The file {sample_fn_curr} should exist"),
          "i" = "Try re running the pipeline from scratch or report the error to GCIMS authors"
        )
      )
    }
    gcimssample <- readRDS(sample_fn_curr)
    # Check for a sample name change:
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
  if (needs_re_saving || is.null(sample_fn_curr)) {
    saveRDS(gcimssample, sample_fn_next)
  } else {
    if (!file.copy(
      sample_fn_curr,
      sample_fn_next
    )) {
      rlang::abort(glue("Could not copy {sample_fn_curr} to {sample_fn_next}."))
    }
  }
  out
}


#' Runs all delayed operations on the object
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @param keep_intermediate A logical, whether to keep the intermediate files of
#' the previous realization once this one finishes. If `NA`, keeping will depend
#' on the `object`.
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
realize <- function(object, keep_intermediate = NA) {
  if (!hasDelayedOps(object)) {
    return(object)
  }

  if (is.na(keep_intermediate)) {
    keep_intermediate <- object@envir$keep_intermediate
  }


  current_dir <- CurrentHashedDir(object)
  next_dir <- NextHashedDir(object)
  dir.create(next_dir, showWarnings = FALSE, recursive = TRUE)

  sample_names <- sampleNames(object)
  names(sample_names) <- sample_names

  pdata <- Biobase::pData(object)
  orig_filenames <- pdata$FileName
  names(orig_filenames) <- sample_names

  extracted_results <- BiocParallel::bplapply(
    X = sample_names,
    FUN = realize_one_sample,
    curr_dir = current_dir,
    next_dir = next_dir,
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
  CurrentHashedDir(object) <- next_dir
  if (!keep_intermediate && !is.null(current_dir)) {
    unlink(current_dir, recursive = TRUE)
  }
  invisible(object)
}
