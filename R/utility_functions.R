require_pkgs <- function(pkg, msgs = NULL, ...) {
    have_pkgs <- purrr::map_lgl(pkg, function(p) {
        requireNamespace(p, quietly = TRUE)
    })
    names(have_pkgs) <- pkg
    if (!all(have_pkgs)) {
        missing_pkgs <- names(have_pkgs)[!have_pkgs]
        parent_call <- format(rlang::caller_call())
        rlang::abort(
            message = c(
                glue::glue("{parent_call} requires additional packages. Please install them. You may want to use:", parent_call = parent_call),
                glue::glue("    BiocManager::install({deparse(missing_pkgs)})", missing_pkgs = missing_pkgs),
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
