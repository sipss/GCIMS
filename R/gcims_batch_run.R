#' Apply a workflow step to all samples
#'
#' @param dir_in Input directory where all .rds files to be processed are
#' @param dir_out Output directory where all .rds filess will be saved
#' @param what Function to run, that takes one sample object and returns one sample object
#' @param ... Arguments to what
#' @param .batch_samples (deprecated) Sample numbers, as expected by other functions in the package, to build `paste0("M", i, ".rds")` functions.
#' @param .batch_returns A function that takes a sample object. `NULL`, by default.
#'
#' @return A list as long as the number of .rds files in `dir_in`. Each list element
#'  is the output of `.batch_returns(obj)`, where obj is the corresponding rds sample
#'  (`NULL` by default).
#' @export
#'
#' @examples
#' # We want to run the following function on one sample, and explore how
#' # the parameter `n` affects the output. Once we feel we have a good `n`,
#' # we want to apply the function to all samples in parallel or whatever.
#'
#' # We have an input directory with two samples
#' dir_in <- tempfile("dir_in")
#' dir.create(dir_in, recursive = TRUE)
#' file.copy(
#'   system.file("extdata", "M3.rds", package = "GCIMS"),
#'   dir_in
#' )
#' file.copy(
#'   system.file("extdata", "M7.rds", package = "GCIMS"),
#'   dir_in
#' )
#'
#' # We have a function we want to try on one sample first, and apply to
#' # all samples afterwards:
#'
#' adds_n <- function(x, n) {
#'   x$data$data_df <- x$data$data_df + n
#'   return(x)
#' }
#'
#' # Run on one sample:
#' one_file <- file.path(dir_in, "M3.rds")
#' one_sample <- readRDS(one_file)
#' # Before:
#' print(one_sample$data$data_df[1:3,1:3])
#' # Compute just on this sample
#' one_sample <- adds_n(one_sample, n  = 4)
#' # After:
#' print(one_sample$data$data_df[1:3,1:3])
#'
#' # We are satisfied with n=4. Apply to all samples:
#' dir_out <- tempfile("dir_out")
#' gcims_batch_run(dir_in, dir_out, adds_n, n=4)
#'
#'
gcims_batch_run <- function(dir_in, dir_out, what, ...,
                            .batch_samples = NULL, .batch_returns = function(x) {NULL}) {
  if (is.null(.batch_samples)) {
    files <- list.files(path = dir_in, pattern = "\\.rds$", recursive = FALSE)
  } else {
    files <- paste0("M", .batch_samples, ".rds")
  }
  if (is.null(names(files))) {
    files <- stats::setNames(files, basename(files))
  }
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  BiocParallel::bplapply(
    files,
    function(f, dir_in, dir_out, what, ...) {
      obj <- readRDS(file.path(dir_in, f))
      obj <- what(obj, ...)
      saveRDS(obj, file.path(dir_out, f))
      .batch_returns(obj)
    },
    dir_in = dir_in,
    dir_out = dir_out,
    what = what,
    ...
  )
}
