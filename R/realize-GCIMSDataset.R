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
    cli_abort("sample_obj should be a GCIMSSample or a file name")
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
      cli_abort(
        message = c(
          "UnexpectedError",
          "x" = "The file {current_filename} should exist",
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
      cli_abort("Could not copy {current_filename} to {next_filename}.")
    }
  }
  res$extracted_objects
}


move_delayed_to_history <- function(object) {
  # Make delayed_ops history and remove them from pending
  object$previous_ops <- c(object$previous_ops, object$delayed_ops)
  object$delayed_ops <- list()
  object
}

aggregate_all_results <- function(object, extracted_results) {
  # Apply to the dataset object
  for (i in seq_along(object$delayed_ops)) {
    delayed_op <- object$delayed_ops[[i]]
    # Extract i-th result for all samples:
    extracted_result <- purrr::map(extracted_results, i)
    # Let the delayed operation aggregate the extracted results and save them in object
    object <- aggregate_result(delayed_op, extracted_result, object)
  }
  object
}

realize_ram <- function(object) {
  sample_objs <- object$samples
  if (is.null(sample_objs)) {
    pdata <- Biobase::pData(object)
    sample_objs <- pdata$FileName
  }
  delayed_ops <- object$delayed_ops
  results <- mymapply(
    FUN = realize_one_sample_ram,
    sample_name = sampleNames(object),
    sample_obj = sample_objs,
    MoreArgs = list(
      delayed_ops = delayed_ops
    ),
    SIMPLIFY = FALSE
  )

  object$samples <- purrr::map(results, "sample_obj")
  extracted_results <- purrr::map(results, "extracted_objects")

  object <- aggregate_all_results(object, extracted_results)
  object <- move_delayed_to_history(object)
  invisible(object)
}

realize_disk <- function(object, keep_intermediate) {
  delayed_ops <- object$delayed_ops
  current_dir <- object$currentHashedDir()
  if (is.null(current_dir)) {
    if (name(delayed_ops[[1]]) != "read_sample") {
      cli_abort(
        message = c(
          "UnexpectedError",
          "x" = "The first operation should have been named read_sample instead of {name(delayed_ops[[1]])}",
          "i" = "This is an unexpected problem. You can try deleting the scratch directory and restart again"
        )
      )
    }
  }
  next_dir <- object$nextHashedDir()
  if (is.null(next_dir)) {
    cli_abort(message = c("UnexpectedError", "x" = "next_dir should not have been NULL"))
  }
  dir.create(next_dir, showWarnings = FALSE, recursive = TRUE)

  sample_names <- sampleNames(object)

  pdata <- Biobase::pData(object)

  if (is.null(current_dir)) {
    current_filenames <- list(NULL)
  } else {
    current_filenames <- file.path(current_dir, paste0(sample_names, ".rds"))
  }

  next_filenames <- file.path(next_dir, paste0(sample_names, ".rds"))

  extracted_results <- mymapply(
    FUN = realize_one_sample_disk,
    sample_name = sample_names,
    orig_filename = pdata$FileName,
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

  object$hasheddir <- basename(next_dir)
  if (is.na(keep_intermediate)) {
    keep_intermediate <- object$keep_intermediate
  }
  if (!keep_intermediate && !is.null(current_dir) && current_dir != next_dir) {
    unlink(current_dir, recursive = TRUE)
  }
  saveRDS(object, file = file.path(next_dir, "GCIMSDataset.rds"))
  object
}


canRealize <- function(object) {
  object$can_realize
}

"canRealize<-" <- function(object, value) {
  value <- as.logical(value)
  if (length(value) != 1 || is.na(value)) {
    cli_abort("canRealize can't be NA or NULL and must be of length one")
  }
  object$can_realize <- value
  object
}

