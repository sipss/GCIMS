setMethod(
  "show",
  "GCIMSSample",
  function(object) {
    axes <- list(
      "drift time" = list(value = dtime(object), unit = "ms"),
      "retention time" = list(value = rtime(object), unit = "s")
    )
    outstring <- "A GCIMS Sample"
    for (axis_name in  names(axes)) {
      axis <- axes[[axis_name]][["value"]]
      axis_unit <- axes[[axis_name]][["unit"]]
      if (length(axis) == 0) {
        first <- NaN
        last <- NaN
        res <- NaN
      } else if (length(axis) == 1) {
        first <- axis[1]
        last <- axis[1]
        res <- NaN
      } else {
        first <- axis[1]
        last <- axis[length(axis)]
        res <- axis[2] - axis[1]
      }
      outstring <- c(
        outstring,
        paste0(" with ", axis_name, " from ", first, " to ", last, " ", axis_unit,
               " (", "step: ", signif(res, digits = 3), " ", axis_unit , ", ",
               "points: ", length(axis), ")")
      )
    }
    # if (length(object@history) > 0) {
    #   outstring <- c(outstring, "History:", sprintf("- %s", object@history))
    # } else {
    #   outstring <- c(outstring, "History:", "No history available")
    # }
    # # FIXME: Give more details about the sample
    cat(paste0(outstring, collapse = "\n"))
  }
)
