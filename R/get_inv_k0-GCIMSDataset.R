methods::setMethod(
  "get_inv_k0",
  "GCIMSDataset",
  function(object){
    delayed_op <- DelayedOperation(
      name = "get_inverse_reduced_mobility",
      fun = get_inv_k0
    )
    object$appendDelayedOp(delayed_op)
    invisible(object)
  }
)
