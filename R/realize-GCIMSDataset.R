realize_one_sample_ram <- function(sample_name, sample_obj, delayed_ops) {
  out <- vector("list", length = length(delayed_ops))
  needs_re_saving <- FALSE
  if (is(sample_obj, "GCIMSSample")) {
    if (!identical(description(sample_obj), sample_name)) {
      # sampleNames were probably updated, files renamed.
      description(sample_obj) <- sample_name
      needs_re_saving <- TRUE
    }
  } else if (is.character(sample_obj)) {
  } else {
    abort("sample_obj should be a GCIMSSample or a file name")
  }
  # Execute:
  for (i in seq_along(delayed_ops)) {
    result <- apply_op_to_sample(delayed_ops[[i]], sample_obj)
    sample_obj <- result$sample
    if (!identical(description(sample_obj), sample_name)) {
      # sampleNames were probably updated, files renamed.
      description(sample_obj) <- sample_name
      needs_re_saving <- TRUE
    }
    if (!is.null(result$extracted_obj)) {
      out[[i]] <- result$extracted_obj
    }
    if (result$needs_resaving) {
      needs_re_saving <- TRUE
    }
  }
  list(
    extracted_objects = out,
    sample_obj = sample_obj,
    needs_re_saving = needs_re_saving
  )
}


realize_one_sample_disk <- function(sample_name, orig_filename, current_filename, next_filename, delayed_ops) {
  if (is.null(current_filename)) {
    # we are going to read_sample, we are starting
    sample_obj <- orig_filename
  } else {
    # We load the file
    if (!file.exists(current_filename)) {
      abort(
        message = c(
          "UnexpectedError",
          "x" = glue("The file {current_filename} should exist"),
          "i" = "Try re running the pipeline from scratch or report the error to GCIMS authors"
        )
      )
    }
    sample_obj <- readRDS(current_filename)
  }
  res <- realize_one_sample_ram(sample_name, sample_obj, delayed_ops)

  # saveRDS, or just copy
  if (res$needs_re_saving || is.null(current_filename)) {
    saveRDS(res$sample_obj, next_filename)
  } else {
    if (!file.copy(
      current_filename,
      next_filename
    )) {
      abort(glue("Could not copy {current_filename} to {next_filename}."))
    }
  }
  res$extracted_objects
}

optimize_delayed_operations <- function(object) {
  if (!hasDelayedOps(object)) {
    return(object)
  }
  delayed_ops <- object@envir$delayed_ops
  # Extra operations that extract the rtime and dtime or TIS and RIC from the object can be delayed
  # Find them:
  where_extract_times <- purrr::map_lgl(delayed_ops, function(op) {name(op) == "extract_dtime_rtime"})
  where_extract_RIC_TIS <- purrr::map_lgl(delayed_ops, function(op) {name(op) == "extract_RIC_and_TIS"})
  where_extract <- where_extract_times | where_extract_RIC_TIS
  # Not found, return:
  if (!any(where_extract)) {
    return(object)
  }
  # Remove those ops:
  delayed_ops[where_extract] <- NULL
  object@envir$delayed_ops <- delayed_ops
  if (any(where_extract_times)) {
    object <- extract_dtime_rtime(object)
  }
  if (any(where_extract_RIC_TIS)) {
    object <- extract_RIC_and_TIS(object)
  }
  object
}

move_delayed_to_history <- function(object) {
  # Make delayed_ops history and remove them from pending
  object@envir$previous_ops <- c(object@envir$previous_ops, object@envir$delayed_ops)
  object@envir$delayed_ops <- list()
  object
}

aggregate_all_results <- function(object, extracted_results) {
  # Apply to the dataset object
  for (i in seq_along(object@envir$delayed_ops)) {
    delayed_op <- object@envir$delayed_ops[[i]]
    # Extract i-th result for all samples:
    extracted_result <- purrr::map(extracted_results, i)
    # Let the delayed operation aggregate the extracted results and save them in object
    object <- aggregate_result(delayed_op, extracted_result, object)
  }
  object
}

runs_in_serial <- function() {
  tryCatch({
    methods::is(BiocParallel::bpparam(), "SerialParam")
  },
  error = function(e) TRUE
  )
}

realize_ram <- function(object) {
  sample_objs <- object@envir$samples
  if (is.null(sample_objs)) {
    pdata <- Biobase::pData(object)
    sample_objs <- file.path(object@envir$base_dir, pdata$FileName)
  }
  delayed_ops <- rlang::as_list(object@envir$delayed_ops)
  if (runs_in_serial()) {
    # I prefer to avoid BiocParallel if running in serial
    mapply_fun <- mapply
  } else {
    mapply_fun <- BiocParallel::bpmapply
  }
  results <- mapply_fun(
    FUN = realize_one_sample_ram,
    sample_name = sampleNames(object),
    sample_obj = sample_objs,
    MoreArgs = list(
      delayed_ops = delayed_ops
    ),
    SIMPLIFY = FALSE
  )

  object@envir$samples <- purrr::map(results, "sample_obj")
  extracted_results <- purrr::map(results, "extracted_objects")

  object <- aggregate_all_results(object, extracted_results)
  object <- move_delayed_to_history(object)
  invisible(object)
}

realize_disk <- function(object, keep_intermediate) {
  delayed_ops <- object@envir$delayed_ops
  current_dir <- CurrentHashedDir(object)
  if (is.null(current_dir)) {
    if (name(delayed_ops[[1]]) != "read_sample") {
      abort(
        message = c(
          "UnexpectedError",
          "x" = glue("The first operation should have been named read_sample instead of {name(delayed_ops[[1]])}"),
          "i" = "This is an unexpected problem. You can try deleting the scratch directory and restart again"
        )
      )
    }
  }
  next_dir <- NextHashedDir(object)
  if (is.null(next_dir)) {
    abort(message = c("UnexpectedError", "x" = "next_dir should not have been NULL"))
  }
  dir.create(next_dir, showWarnings = FALSE, recursive = TRUE)

  sample_names <- sampleNames(object)

  pdata <- Biobase::pData(object)
  orig_filenames <- file.path(object@envir$base_dir, pdata$FileName)

  if (is.null(current_dir)) {
    current_filenames <- list(NULL)
  } else {
    current_filenames <- file.path(current_dir, paste0(sample_names, ".rds"))
  }

  next_filenames <- file.path(next_dir, paste0(sample_names, ".rds"))

  if (runs_in_serial()) {
    # I prefer to avoid BiocParallel if running in serial
    mapply_fun <- mapply
  } else {
    mapply_fun <- BiocParallel::bpmapply
  }

  extracted_results <- mapply_fun(
    FUN = realize_one_sample_disk,
    sample_name = sample_names,
    orig_filename = orig_filenames,
    current_filename = current_filenames,
    next_filename = next_filenames,
    MoreArgs = list(
      delayed_ops = delayed_ops
    ),
    SIMPLIFY = FALSE
  )

  object <- aggregate_all_results(object, extracted_results)
  # Make delayed_ops history and remove them from pending
  object <- move_delayed_to_history(object)

  CurrentHashedDir(object) <- next_dir
  if (is.na(keep_intermediate)) {
    keep_intermediate <- object@envir$keep_intermediate
  }
  if (!keep_intermediate && !is.null(current_dir)) {
    unlink(current_dir, recursive = TRUE)
  }
  saveRDS(object, file = file.path(next_dir, "GCIMSDataset.rds"))
  object
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
setMethod("realize", "GCIMSDataset",  function(object, keep_intermediate = NA) {
  if (!hasDelayedOps(object)) {
    return(object)
  }

  if (!canRealize(object)) {
    abort(
      message = c(
        "UnexpectedError",
        paste0(
          "Tried to realize a dataset that could not be realized. This is ",
          "a bug in the package. Please report it to https://github.com/sipss/GCIMS/issues"
        )
      )
    )
  }
  canRealize(object) <- FALSE
  on.exit({canRealize(object) <- TRUE})

  object <- optimize_delayed_operations(object)
  if (isTRUE(object@envir$on_ram)) {
    object <- realize_ram(object)
  } else {
    object <- realize_disk(object, keep_intermediate = keep_intermediate)
  }
  invisible(object)
})


canRealize <- function(object) {
  object@envir$can_realize
}

"canRealize<-" <- function(object, value) {
  value <- as.logical(value)
  if (length(value) != 1 || is.na(value)) {
    abort("canRealize can't be NA or NULL and must be of length one")
  }
  object@envir$can_realize <- value
  object
}

