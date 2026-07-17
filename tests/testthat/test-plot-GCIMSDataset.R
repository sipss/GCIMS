# --- resolve_intensity_range() -----------------------------------------

make_counting_range <- function(value) {
  n_calls <- 0L
  fn <- function() {
    n_calls <<- n_calls + 1L
    value
  }
  attr(fn, "n_calls") <- function() n_calls
  fn
}

test_that("resolve_intensity_range('global'/'ranged') calls only the matching closure, once", {
  global_range <- make_counting_range(c(0, 10))
  ranged_range <- make_counting_range(c(2, 8))

  expect_equal(resolve_intensity_range("global", global_range, ranged_range), c(0, 10))
  expect_equal(attr(global_range, "n_calls")(), 1L)
  expect_equal(attr(ranged_range, "n_calls")(), 0L)
})

test_that("resolve_intensity_range() with a fixed numeric vector never calls either closure", {
  global_range <- make_counting_range(c(0, 10))
  ranged_range <- make_counting_range(c(2, 8))

  expect_equal(resolve_intensity_range(c(1, 9), global_range, ranged_range), c(1, 9))
  expect_equal(attr(global_range, "n_calls")(), 0L)
  expect_equal(attr(ranged_range, "n_calls")(), 0L)
})

test_that("resolve_intensity_range() resolves mixed list(min=, max=) elements independently", {
  global_range <- make_counting_range(c(0, 10))
  ranged_range <- make_counting_range(c(2, 8))

  expect_equal(resolve_intensity_range(list(min = 1, max = "global"), global_range, ranged_range), c(1, 10))
  expect_equal(attr(global_range, "n_calls")(), 1L)
  expect_equal(attr(ranged_range, "n_calls")(), 0L)

  expect_equal(resolve_intensity_range(list("ranged", 99), global_range, ranged_range), c(2, 99))
  expect_equal(attr(ranged_range, "n_calls")(), 1L)
})

test_that("resolve_intensity_range() rejects invalid specs", {
  global_range <- make_counting_range(c(0, 10))
  ranged_range <- make_counting_range(c(2, 8))

  expect_error(resolve_intensity_range("nope", global_range, ranged_range), "intensity_range should be")
  expect_error(resolve_intensity_range(c(1, 2, 3), global_range, ranged_range), "intensity_range should be")
  expect_error(resolve_intensity_range(list("nope", 1), global_range, ranged_range), "should be a number")
})

# --- resolve_page_grid() --------------------------------------------------

test_that("resolve_page_grid() with neither nrow nor ncol given follows the exact-fit-then-3x3-cap table", {
  expected <- list(
    `1` = c(1, 1), `2` = c(1, 2), `3` = c(1, 3), `4` = c(2, 2),
    `5` = c(2, 3), `6` = c(2, 3), `7` = c(3, 3), `9` = c(3, 3), `25` = c(3, 3)
  )
  for (n in names(expected)) {
    g <- resolve_page_grid(NULL, NULL, as.integer(n))
    expect_equal(c(g$nrow, g$ncol), expected[[n]], info = n)
  }
})

test_that("resolve_page_grid() with exactly one of nrow/ncol given fills the other", {
  # num_samples <= 9: fit everyone on one page
  expect_equal(resolve_page_grid(NULL, 2, 5)$nrow, 3)
  expect_equal(resolve_page_grid(4, NULL, 9)$ncol, 3)

  # num_samples > 9: missing dimension is fixed at 3, regardless of the given one
  expect_equal(resolve_page_grid(NULL, 5, 20)$nrow, 3)
  expect_equal(resolve_page_grid(7, NULL, 50)$ncol, 3)
})

test_that("resolve_page_grid() with both given uses them as-is", {
  g <- resolve_page_grid(4, 5, 100)
  expect_equal(c(g$nrow, g$ncol), c(4, 5))
})

# --- intensity_range caching ----------------------------------------------

make_range_dataset <- function() {
  s1 <- GCIMSSample(drift_time = 1:5, retention_time = 1:5, data = matrix(1:25, nrow = 5))
  s2 <- GCIMSSample(drift_time = 1:5, retention_time = 1:5, data = matrix(100:124, nrow = 5))
  GCIMSDataset$new_from_list(samples = list(s1 = s1, s2 = s2), on_ram = TRUE, scratch_dir = NULL)
}

test_that("realizing a dataset populates intensity_range alongside TIS/RIC", {
  old_bpparam <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(old_bpparam))

  ds <- make_range_dataset()
  ds$realize()

  expect_equal(unname(ds$intensity_range["s1", ]), c(1, 25))
  expect_equal(unname(ds$intensity_range["s2", ]), c(100, 124))
})

test_that("subsetting a dataset resets the cached intensity_range", {
  old_bpparam <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(old_bpparam))

  ds <- make_range_dataset()
  ds$realize()
  expect_false(is.null(ds$intensity_range))

  ds$subset("s1", inplace = TRUE)

  expect_null(ds$intensity_range)
})

# --- plot(GCIMSDataset) ----------------------------------------------------

test_that("plot(GCIMSDataset) returns a combined plot using the cached global range by default", {
  ds <- make_range_dataset()

  p <- plot(ds)

  expect_s3_class(p, "ggplot")
})

test_that("plot(GCIMSDataset, sample =) restricts which samples are plotted", {
  ds <- make_range_dataset()

  expect_no_error(plot(ds, sample = "s1"))
  expect_no_error(plot(ds, sample = 1))
})

test_that("plot(GCIMSDataset, intensity_range = 'global') computes a baseline-removed range with remove_baseline = TRUE", {
  ds <- make_range_dataset()
  ds <- estimateBaseline(ds, dt_peak_fwhm_ms = 1, dt_region_multiplier = 2, rt_length_s = 2)
  ds$realize()

  expect_no_error(plot(ds, remove_baseline = TRUE))
})

test_that("plot(GCIMSDataset, intensity_range = 'ranged') works with remove_baseline = TRUE", {
  ds <- make_range_dataset()
  ds <- estimateBaseline(ds, dt_peak_fwhm_ms = 1, dt_region_multiplier = 2, rt_length_s = 2)
  ds$realize()

  expect_no_error(plot(ds, remove_baseline = TRUE, intensity_range = "ranged"))
})

test_that("plot(GCIMSDataset, intensity_range = c(min, max)) accepts fixed limits without touching the dataset", {
  ds <- make_range_dataset()

  expect_no_error(plot(ds, intensity_range = c(0, 200)))
})

# --- pagination -------------------------------------------------------------

make_paged_dataset <- function(n = 11) {
  samples <- stats::setNames(
    purrr::map(seq_len(n), function(i) {
      GCIMSSample(drift_time = 1:5, retention_time = 1:5, data = matrix(rep(i, 25), nrow = 5))
    }),
    paste0("s", seq_len(n))
  )
  GCIMSDataset$new_from_list(samples = samples, on_ram = TRUE, scratch_dir = NULL)
}

test_that("plot(GCIMSDataset) with more than 9 samples defaults to a 3x3 first page", {
  ds <- make_paged_dataset(11)

  expect_no_error(plot(ds))
  expect_no_error(plot(ds, page = 1))
  expect_no_error(plot(ds, page = 2))
})

test_that("plot(GCIMSDataset, page =) errors clearly when out of bounds", {
  ds <- make_paged_dataset(11)

  expect_error(plot(ds, page = 3), "page.*3.*out of bounds")
  expect_error(plot(ds, page = 0), "page.*0.*out of bounds")
})

test_that("plot(GCIMSDataset, intensity_range = 'ranged') is computed over every selected sample, not just the current page", {
  # 11 samples, values 1..11 -> page 1 (samples 1-9) would only see 1..9 if
  # ranged were page-scoped, but sample 11 (value 11) is on page 2.
  ds <- make_paged_dataset(11)

  captured <- list()
  testthat::local_mocked_bindings(
    mat_to_gplot = function(intmat, ..., intensity_range = NULL) {
      captured[[length(captured) + 1]] <<- intensity_range
      ggplot2::ggplot()
    },
    .package = "GCIMS"
  )

  plot(ds, intensity_range = "ranged", page = 1)
  range_page1 <- captured[[1]]
  captured <<- list()
  plot(ds, intensity_range = "ranged", page = 2)
  range_page2 <- captured[[1]]

  expect_equal(range_page1, c(1, 11))
  expect_equal(range_page2, c(1, 11))
})
