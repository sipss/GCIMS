#' Filter GCIMSDataset samples by retention time
#'
#' @param object A [GCIMSDataset] object
#' @inheritParams dt_rt_range_normalization
#' @return The given object, with a delayed operation to filter retention times
#' @export
#' @examples
#' base_dir <- system.file("extdata", "sample_formats", package = "GCIMS")
#' annot <- data.frame(SampleID = "Sample1", FileName = "small.mea.gz")
#' dataset <- GCIMSDataset(annot, base_dir)
#' filterRt(dataset, rt_range = c(5, 50))
setMethod(
  "filterRt",
  "GCIMSDataset",
  function(object, rt_range) {
    delayed_op <- GCIMSDelayedOp(
      name = "filterRt",
      fun = filterRt,
      params = list(rt_range = rt_range)
    )
    object <- appendDelayedOp(object, delayed_op)
    # We recompute these, but  maybe we could just reset them to zero...
    object <- extract_dtime_rtime(object)
    object <- extract_RIC_and_TIS(object)
    invisible(object)
  }
)

#' Filter GCIMSDataset samples by drift time
#'
#' @param object A [GCIMSDataset] object
#' @inheritParams dt_rt_range_normalization
#' @export
#' @return The given object, with a delayed operation to filter retention times
#' @examples
#' base_dir <- system.file("extdata", "sample_formats", package = "GCIMS")
#' annot <- data.frame(SampleID = "Sample1", FileName = "small.mea.gz")
#' dataset <- GCIMSDataset(annot, base_dir)
#' filterDt(dataset, dt_range = c(5, 10))
setMethod(
  "filterDt",
  "GCIMSDataset",
  function(object, dt_range) {
    delayed_op <- GCIMSDelayedOp(
      name = "filterDt",
      fun = filterDt,
      params = list(dt_range = dt_range)
    )
    object <- appendDelayedOp(object, delayed_op)
    # We recompute these, but  maybe we could just reset them to zero...
    object <- extract_dtime_rtime(object)
    object <- extract_RIC_and_TIS(object)
    invisible(object)
  }
)

