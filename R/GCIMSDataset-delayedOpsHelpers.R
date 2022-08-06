hasDelayedOps <- function(object) {
  length(object@envir$delayed_ops) > 0
}

appendDelayedOp <- function(object, op) {
  object@envir$delayed_ops <- c(object@envir$delayed_ops, list(op))
  object
}
