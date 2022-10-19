

new_progress_bar <- function(...) {
  if (!requireNamespace("progress", quietly = TRUE)) {
    inform(
      message = c("i" = 'Use install.packages("progress") to get a progress bar'),
      class = "GCIMS_suggest_install_progress",
      .frequency = "once",
      .frequeny_id = "suggest_progress_installation"
    )
    # This code path returns a dummy object exposing a compatible API of
    # progress:progress_bar$new(...).
    dummy_pbar <- list(
      tick = function() {}
    )
    return(dummy_pbar)
  }
  progress::progress_bar$new(...)
}

require_pkgs <- function(pkg, msgs = NULL, ...) {
  have_pkgs <- purrr::map_lgl(pkg, function(p) {requireNamespace(p, quietly = TRUE)})
  names(have_pkgs) <- pkg
  if (!all(have_pkgs)) {
    missing_pkgs <- names(have_pkgs)[!have_pkgs]
    aval_pkgs <- rownames(utils::available.packages())
    missing_cran_pkgs <- intersect(missing_pkgs, aval_pkgs)
    missing_bioc_pkgs <- setdiff(missing_pkgs, missing_cran_pkgs)
    if (length(missing_bioc_pkgs) > 0) {
      if (!"BiocManager" %in% rownames(utils::installed.packages())) {
        missing_cran_pkgs <- c(missing_cran_pkgs, "BiocManager")
      }
    }
    and_string <- character()
    if (length(missing_cran_pkgs) > 0) {
      missing_cran_pkgs <- deparse(missing_cran_pkgs)
      and_string <- ""
    }
    if (length(missing_bioc_pkgs) > 0) {
      missing_bioc_pkgs <- deparse(missing_bioc_pkgs)
      and_string <- " and"
    }
    parent_call <- format(rlang::caller_call())
    abort(
      message = c(
        glue::glue("{parent_call} requires additional packages. Please install them. You may want to use:", parent_call = parent_call),
        glue::glue("    install.packages({missing_cran_pkgs}){and_string}", missing_cran_pkgs = missing_cran_pkgs, and_string = and_string),
        glue::glue("    BiocManager::install({missing_bioc_pkgs})", missing_bioc_pkgs = missing_bioc_pkgs),
        msgs
      ),
      ...
    )
  }
}


units_to_points <- function(length_phys, step_phys, must_odd = FALSE) {
  length_pts <- round(length_phys/step_phys)
  if (must_odd) {
    length_pts <- length_pts + (length_pts %% 2 == 0) # the filter length in points
  }
  length_pts
}
