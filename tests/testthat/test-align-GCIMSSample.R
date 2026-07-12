bump <- function(rt, center, height, sd) height * exp(-(rt - center)^2 / (2 * sd^2))

# A GCIMSSample with a single informative drift-time row (the RIP row), so
# TIS/RIC-based alignment logic has an unambiguous target to lock onto.
make_rip_sample <- function(dt, rt, rip_idx, rip_baseline = 100, spike_col = NULL, spike_height = 1e5) {
  data <- matrix(1, nrow = length(dt), ncol = length(rt))
  data[rip_idx, ] <- rip_baseline
  if (!is.null(spike_col)) {
    data[rip_idx, spike_col] <- spike_height
  }
  GCIMSSample(drift_time = dt, retention_time = rt, data = data)
}

# A GCIMSSample whose RIP row is `bump`-shaped (inverted, as getRIC() expects),
# so its RIC is a clean, controllable Gaussian bump centered at `center`.
make_bump_sample <- function(dt, rt, center, height = 900, sd = 3) {
  data <- matrix(1, nrow = length(dt), ncol = length(rt))
  data[1, ] <- 1000 - bump(rt, center = center, height = height, sd = sd)
  GCIMSSample(drift_time = dt, retention_time = rt, data = data)
}

test_that("alignDt computes the multiplicative correction from the RIP position", {
  dt <- seq(1, 10, by = 0.5)
  rt <- 1:20
  rip_idx <- 10 # dt[10] == 5.5
  s <- make_rip_sample(dt, rt, rip_idx)
  expect_equal(dt[which.max(getTIS(s))], dt[rip_idx])

  rip_ref_ms <- 6
  aligned <- alignDt(s, rip_ref_ms = rip_ref_ms)

  expect_equal(aligned@proc_params$align$dt_kcorr, rip_ref_ms / dt[rip_idx])
  expect_equal(dim(intensity(aligned)), dim(intensity(s)))
})

test_that("alignDt with rip_ref_ms already matching the RIP position is a no-op correction", {
  dt <- seq(1, 10, by = 0.5)
  rt <- 1:20
  rip_idx <- 10
  s <- make_rip_sample(dt, rt, rip_idx)

  aligned <- alignDt(s, rip_ref_ms = dt[rip_idx])

  expect_equal(aligned@proc_params$align$dt_kcorr, 1)
})

test_that("alignRt_ip slices the sample to rt_ref, starting at the injection point", {
  dt <- 1:3
  rt <- 1:30
  injection_col <- 10
  s <- make_rip_sample(dt, rt, rip_idx = 2, spike_col = injection_col)
  expect_equal(unname(which.min(getRIC(s))), injection_col)

  min_start <- 3
  rt_ref <- seq(100, 114) # length 15
  aligned <- alignRt_ip(s, min_start = min_start, rt_ref = rt_ref)

  expect_identical(rtime(aligned), rt_ref)
  expected_start <- injection_col - min_start
  expect_equal(
    unname(intensity(aligned)),
    unname(intensity(s))[, expected_start:(expected_start + length(rt_ref) - 1)]
  )
})

test_that("prealign combines drift-time and injection-point alignment without introducing NAs", {
  dt <- seq(1, 10, by = 0.5)
  rt <- 1:30
  rip_idx <- 10
  injection_col <- 10
  s <- make_rip_sample(dt, rt, rip_idx = rip_idx, spike_col = injection_col)

  rip_ref_ms <- 6
  min_start <- 3
  rt_ref <- seq(100, 114)

  aligned <- prealign(
    s,
    align_dt = TRUE, align_ip = TRUE,
    rip_ref_ms = rip_ref_ms, min_start = min_start, rt_ref = rt_ref
  )

  expect_equal(aligned@proc_params$align$dt_kcorr, rip_ref_ms / dt[rip_idx])
  expect_identical(rtime(aligned), rt_ref)
  expect_false(anyNA(intensity(aligned)))
})

test_that("prealign with align_dt = FALSE, align_ip = FALSE leaves axes untouched", {
  dt <- seq(1, 10, by = 0.5)
  rt <- 1:30
  s <- make_rip_sample(dt, rt, rip_idx = 10, spike_col = 10)

  aligned <- prealign(
    s,
    align_dt = FALSE, align_ip = FALSE,
    rip_ref_ms = 6, min_start = 3, rt_ref = seq(100, 114)
  )

  expect_equal(aligned@proc_params$align$dt_kcorr, 1)
  expect_identical(dtime(aligned), dt)
  expect_identical(rtime(aligned), rt)
})

test_that("prealign refuses to align a sample that is already entirely missing", {
  dt <- 1:2
  rt <- 1:10
  s_na <- GCIMSSample(drift_time = dt, retention_time = rt, data = matrix(NA_real_, nrow = 2, ncol = 10))

  expect_error(
    prealign(s_na, align_dt = TRUE, align_ip = TRUE, rip_ref_ms = 1, min_start = 1, rt_ref = 1:5),
    "Align is impossible"
  )
})

test_that("align(method_rt = 'ptw') interpolates the reference RIC when it's on a different rt grid", {
  dt <- 1:2
  rt <- 1:60

  ref_sample <- make_bump_sample(dt, rt, center = 30)
  ric_ref <- getRIC(ref_sample)
  ric_ref_rt_diff <- rt + 0.5 # different grid than the sample's own rt, still overlapping

  shifted_sample <- make_bump_sample(dt, rt, center = 24)
  aligned <- align(shifted_sample, method_rt = "ptw", ric_ref = ric_ref, ric_ref_rt = ric_ref_rt_diff)

  expect_false(anyNA(intensity(aligned)))
})

test_that("alignRt_ptw leaves an already well-aligned sample uncorrected (polynomial order 0)", {
  dt <- 1:2
  rt <- 1:60

  ref_sample <- make_bump_sample(dt, rt, center = 30)
  ric_ref <- getRIC(ref_sample)
  same_sample <- make_bump_sample(dt, rt, center = 30)

  aligned <- alignRt_ptw(same_sample, ric_ref = ric_ref, ric_ref_rt = rt)

  expect_identical(aligned@proc_params$align$w, seq_along(rt))
})

test_that("align(method_rt = 'ptw') warps a shifted retention-time profile onto the reference", {
  dt <- 1:2
  rt <- 1:60

  ref_sample <- make_bump_sample(dt, rt, center = 30)
  ric_ref <- getRIC(ref_sample)

  shifted_sample <- make_bump_sample(dt, rt, center = 24)
  aligned <- align(shifted_sample, method_rt = "ptw", ric_ref = ric_ref, ric_ref_rt = rt)

  ric_after <- getRIC(aligned)
  expect_false(anyNA(intensity(aligned)))
  expect_equal(rtime(aligned)[which.max(ric_after)], 30)

  ric_ref_interp <- signal::interp1(rt, ric_ref, xi = rtime(aligned), extrap = FALSE)
  expect_equal(cor(ric_ref_interp, ric_after, use = "complete.obs"), 1, tolerance = 1e-6)
})

test_that("align(method_rt = 'pow') fails with a clear message when pow isn't installed", {
  skip_if(requireNamespace("pow", quietly = TRUE), "pow is installed; not testing the missing-dependency guard")

  dt <- 1:2
  rt <- 1:60
  s <- make_bump_sample(dt, rt, center = 30)
  ric_ref <- getRIC(s)

  expect_error(
    align(s, method_rt = "pow", ric_ref = ric_ref, ric_ref_rt = rt),
    "pow.*is not on CRAN/Bioconductor"
  )
})

test_that("align aborts instead of silently returning an empty/all-NA sample", {
  dt <- 1:2
  rt <- 1:60
  s <- make_bump_sample(dt, rt, center = 30)
  ric_ref <- getRIC(s)
  disjoint_ric_ref_rt <- rt + 10000 # no overlap with the sample's own retention time

  expect_error(
    align(s, method_rt = "none", ric_ref = ric_ref, ric_ref_rt = disjoint_ric_ref_rt),
    "missing values"
  )
})

test_that("align aborts when the data is already entirely missing before subsetting", {
  s <- GCIMSSample(drift_time = 1:2, retention_time = 1:10, data = matrix(NA_real_, nrow = 2, ncol = 10))

  expect_error(
    align(s, method_rt = "none", ric_ref = rep(1, 10), ric_ref_rt = 1:10),
    "After aligning retention time"
  )
})

test_that("alignDt truncates dt_max_ms when the correction shrinks the high end of the drift time axis", {
  dt <- seq(1, 10, by = 0.5)
  rt <- 1:20
  rip_idx <- 10 # dt[10] == 5.5
  s <- make_rip_sample(dt, rt, rip_idx)

  # rip_ref_ms < dt[rip_idx] gives Kcorr < 1, so dt_corr shrinks below dt_end:
  aligned <- alignDt(s, rip_ref_ms = 3)

  expect_lt(aligned@proc_params$align$dt_kcorr, 1)
  expect_lt(aligned@proc_params$align$dt_max_ms, max(dt))
})

test_that("prealign aborts when drift-time alignment yields an all-NA sample", {
  # alignDt() cannot naturally produce an all-NA sample through the public
  # API (it interpolates with extrap = TRUE), so this defensive check is
  # forced by mocking alignDt() directly.
  testthat::local_mocked_bindings(
    alignDt = function(object, rip_ref_ms) {
      object@data[] <- NA_real_
      object@proc_params$align <- list(
        dt_kcorr = 1, dt_min_ms = min(dtime(object)), dt_max_ms = max(dtime(object))
      )
      object
    },
    .package = "GCIMS"
  )
  s <- make_rip_sample(1:2, 1:10, rip_idx = 1)

  expect_error(
    prealign(s, align_dt = TRUE, align_ip = FALSE, rip_ref_ms = 1, min_start = 1, rt_ref = 1:5),
    "After aligning drift times"
  )
})

test_that("prealign aborts when the post-alignment subset leaves no data", {
  # None of prealign()'s own dt_range/rt_range computations can naturally
  # produce an empty/all-NA subset through the public API, so this
  # defensive check is forced by mocking subset.GCIMSSample() directly.
  testthat::local_mocked_bindings(
    subset.GCIMSSample = function(x, ...) {
      x@data[] <- NA_real_
      x
    },
    .package = "GCIMS"
  )
  s <- make_rip_sample(1:2, 1:10, rip_idx = 1)

  expect_error(
    prealign(s, align_dt = FALSE, align_ip = FALSE, rip_ref_ms = 1, min_start = 1, rt_ref = 1:5),
    "After aligning and subsetting"
  )
})

test_that("alignRt_ptw truncates rt_max_s when the correction shrinks the high end of the retention time axis", {
  dt <- 1:2
  rt <- 1:60

  ref_sample <- make_bump_sample(dt, rt, center = 30)
  ric_ref <- getRIC(ref_sample)
  shifted_sample <- make_bump_sample(dt, rt, center = 45)

  aligned <- alignRt_ptw(shifted_sample, ric_ref = ric_ref, ric_ref_rt = rt)

  expect_lt(aligned@proc_params$align$rt_max_s, max(rt))
})
