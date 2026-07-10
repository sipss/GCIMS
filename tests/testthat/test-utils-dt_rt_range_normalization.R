dt <- seq(1, 10, by = 1) # 1..10, step 1
rt <- seq(0, 5, by = 0.5) # 0..5, step 0.5

test_that("dt_rt_range_normalization with a length-2 range selects the inclusive interval", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt, dt_range = c(3, 6))

  expect_equal(out$dt_ms_min, 3)
  expect_equal(out$dt_ms_max, 6)
  expect_equal(out$dt_idx_min, 3L)
  expect_equal(out$dt_idx_max, 6L)
  expect_equal(out$dt_logical, dt >= 3 & dt <= 6)
})

test_that("dt_rt_range_normalization with a single value selects the exact match", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt, dt_range = 5)

  expect_equal(out$dt_ms_min, 5)
  expect_equal(out$dt_ms_max, 5)
  expect_equal(out$dt_idx_min, 5L)
  expect_equal(out$dt_idx_max, 5L)
  expect_equal(which(out$dt_logical), 5L)
})

test_that("dt_rt_range_normalization with a single value accepts one time step of rounding tolerance", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt, dt_range = 5.5) # dt step is 1

  expect_equal(out$dt_idx_min, 5L) # closest point
})

test_that("dt_rt_range_normalization with a single value out of range errors", {
  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, dt_range = 50),
    "not in"
  )
})

test_that("dt_rt_range_normalization with a logical dt_idx uses it directly", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt, dt_idx = dt >= 3 & dt <= 6)

  expect_equal(out$dt_idx_min, 3L)
  expect_equal(out$dt_idx_max, 6L)
  expect_equal(out$dt_ms_min, 3)
  expect_equal(out$dt_ms_max, 6)
})

test_that("dt_rt_range_normalization with a logical dt_idx of the wrong length errors", {
  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, dt_idx = c(TRUE, FALSE)),
    "length"
  )
})

test_that("dt_rt_range_normalization with a numeric dt_idx accepts non-contiguous indices", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt, dt_idx = c(2, 4, 8))

  expect_equal(out$dt_idx_min, 2)
  expect_equal(out$dt_idx_max, 8)
  expect_equal(which(out$dt_logical), c(2L, 4L, 8L))
})

test_that("dt_rt_range_normalization with an out-of-range numeric dt_idx errors", {
  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, dt_idx = c(1, 100)),
    "out of range"
  )
})

test_that("dt_rt_range_normalization with neither dt_range nor dt_idx defaults to the full axis", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt)

  expect_equal(out$dt_idx_min, 1L)
  expect_equal(out$dt_idx_max, length(dt))
  expect_true(all(out$dt_logical))
})

test_that("dt_rt_range_normalization errors when both dt_range and dt_idx are given", {
  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, dt_range = c(1, 2), dt_idx = c(1, 2)),
    "either dt_range or dt_idx"
  )
})

test_that("dt_rt_range_normalization reuses an already-normalized dt_range as-is", {
  already <- dt_rt_range_normalization(dt = dt, rt = rt, dt_range = c(3, 6))

  out <- dt_rt_range_normalization(dt = dt, rt = rt, dt_range = already)

  expect_identical(out, already)
})

test_that("dt_rt_range_normalization with an empty axis and no range/idx returns empty results", {
  out <- dt_rt_range_normalization(dt = numeric(0), rt = rt)

  expect_equal(out$dt_idx_min, integer(0))
  expect_equal(out$dt_logical, logical(0))
})

test_that("dt_rt_range_normalization rejects a dt_range of invalid length", {
  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, dt_range = c(1, 2, 3)),
    "length 2 or a single number"
  )
})

test_that("dt_rt_range_normalization applies the same logic symmetrically to rt", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt, rt_range = c(1, 2))

  expect_equal(out$rt_s_min, 1)
  expect_equal(out$rt_s_max, 2)
  expect_equal(out$rt_logical, rt >= 1 & rt <= 2)

  out_idx <- dt_rt_range_normalization(dt = dt, rt = rt, rt_idx = c(2, 4))
  expect_equal(out_idx$rt_idx_min, 2)
  expect_equal(out_idx$rt_idx_max, 4)

  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, rt_range = c(1, 2), rt_idx = c(1, 2)),
    "either rt_range or rt_idx"
  )
})

test_that("dt_rt_range_normalization with a single rt value selects the closest match within tolerance", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt, rt_range = 2) # exact match, rt step is 0.5

  expect_equal(out$rt_s_min, 2)
  expect_equal(out$rt_idx_min, which(rt == 2))

  out_tol <- dt_rt_range_normalization(dt = dt, rt = rt, rt_range = 2.2)
  expect_equal(out_tol$rt_idx_min, which(rt == 2))

  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, rt_range = 500),
    "not in"
  )
})

test_that("dt_rt_range_normalization with a logical rt_idx validates its length", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt, rt_idx = rt >= 1 & rt <= 2)

  expect_equal(out$rt_s_min, 1)
  expect_equal(out$rt_s_max, 2)

  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, rt_idx = c(TRUE, FALSE)),
    "length"
  )
})

test_that("dt_rt_range_normalization rejects an out-of-range numeric rt_idx", {
  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, rt_idx = c(1, 1000)),
    "out of range"
  )
})

test_that("dt_rt_range_normalization defaults rt to the full axis, or empty for an empty axis", {
  out <- dt_rt_range_normalization(dt = dt, rt = rt)
  expect_equal(out$rt_idx_min, 1L)
  expect_equal(out$rt_idx_max, length(rt))
  expect_true(all(out$rt_logical))

  out_empty <- dt_rt_range_normalization(dt = dt, rt = numeric(0))
  expect_equal(out_empty$rt_idx_min, integer(0))
  expect_equal(out_empty$rt_logical, logical(0))
})

test_that("dt_rt_range_normalization rejects an rt_range of invalid length", {
  expect_error(
    dt_rt_range_normalization(dt = dt, rt = rt, rt_range = c(1, 2, 3)),
    "length 2 or a single number"
  )
})
