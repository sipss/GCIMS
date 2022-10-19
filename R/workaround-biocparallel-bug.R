# This whole file is a workaround for:
# - https://github.com/Bioconductor/BiocParallel/pull/227
#
# Please once that's merged and released, ensure you have
# BiocParallel (>= 1.xx.xx?)
# in the DESCRIPTION file, replace all usages of mybpmapply() with
# BiocParallel::bpmapply() and delete this file.
#
# Some of these functions are copies from the BiocParallel package
.getDotsForMapply <- function (...)
{
  ddd <- list(...)
  if (!length(ddd))
    return(list(list()))
  len <- vapply(ddd, length, integer(1L))
  if (!all(len == len[1L])) {
    max.len <- max(len)
    if (max.len && any(len == 0L))
      stop("zero-length and non-zero length inputs cannot be mixed")
    if (any(max.len%%len))
      warning("longer argument not a multiple of length of vector")
    ddd <- lapply(ddd, rep_len, length.out = max.len)
  }
  ddd
}

.mrename <- function (results, dots, USE.NAMES = TRUE)
{
  if (USE.NAMES) {
    if (length(dots))
      dots <- dots[[1L]]
    if (is.character(dots) && is.null(names(dots))) {
      names(results) <- dots
    }
    else {
      names(results) <- names(dots)
    }
  }
  else {
    results <- unname(results)
  }
  results
}

.simplify <- function (results, SIMPLIFY = FALSE)
{
  if (SIMPLIFY && length(results))
    results <- simplify2array(results)
  results
}

.wrap <- function(.i, .FUN, .ddd, .MoreArgs) {
  dots <- lapply(.ddd, `[`, .i)
  .mapply(.FUN, dots, .MoreArgs)[[1L]]
}


mybpmapply <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE,
                     BPREDO = list(), BPPARAM = BiocParallel::bpparam(), BPOPTIONS = BiocParallel::bpoptions()) {
  ddd <- .getDotsForMapply(...)
  if (!length(ddd) || !length(ddd[[1L]]))
    return(.mrename(list(), ddd, USE.NAMES))

  FUN <- match.fun(FUN)

  res <- BiocParallel:::bplapply(
    X = seq_along(ddd[[1L]]),
    .wrap,
    .FUN = FUN,
    .ddd = ddd,
    .MoreArgs = MoreArgs,
    BPREDO=BPREDO,
    BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS
  )
  .simplify(.mrename(res, ddd, USE.NAMES), SIMPLIFY)
}
