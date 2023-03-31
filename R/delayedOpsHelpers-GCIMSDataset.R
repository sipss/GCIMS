hasDelayedOps <- function(object) {
  cli_warn(c("Deprecated. Use {.code object$hasDelayedOps()} instead."))
  object$hasDelayedOps()
}

appendDelayedOp <- function(object, op) {
  cli_warn(c("Deprecated. Use {.code object$appendDelayedOp(op)} instead."))
  object$appendDelayedOp(op)
  object
}
