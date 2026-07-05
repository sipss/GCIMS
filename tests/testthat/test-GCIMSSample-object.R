test_that("Can create minimal GCIMSSample object", {
  # dummy values:
  dt <- 1:2
  rt <- 1:3
  data <- matrix(seq(length(dt)*length(rt)), nrow= length(dt), ncol = length(rt))
  obj <- GCIMSSample(drift_time = dt, retention_time = rt, data=data)
  expect_s4_class(obj, "GCIMSSample")
  expect_equal(slot(obj, "class_version"), .CURRENT_GCIMSSAMPLE_CLASS_VERSION)
  expect_equal(dtime(obj), dt)
  expect_equal(rtime(obj), rt)
  expect_equal(unname(intensity(obj)), data)
})


test_that("Can create GCIMSSample object with optional slot", {
  # dummy values:
  dt <- 1:2
  rt <- 1:3
  data <- matrix(seq(length(dt)*length(rt)), nrow= length(dt), ncol = length(rt))
  drift_gas <- "nitrogen"
  obj <- GCIMSSample(drift_time = dt, retention_time = rt, data=data, drift_gas = drift_gas)
  expect_s4_class(obj, "GCIMSSample")
  expect_equal(slot(obj, "class_version"), .CURRENT_GCIMSSAMPLE_CLASS_VERSION)
  expect_equal(dtime(obj), dt)
  expect_equal(rtime(obj), rt)
  expect_equal(unname(intensity(obj)), data)
  expect_equal(slot(obj, "drift_gas"), drift_gas)
})

test_that("Invalid argument gives an error", {
  # dummy values:
  dt <- 1:2
  rt <- 1:3
  data <- matrix(seq(length(dt)*length(rt)), nrow= length(dt), ncol = length(rt))
  drift_gas <- "nitrogen"
  expect_error(
    GCIMSSample(drift_time = dt, retention_time = rt, data=data, wrong_argument = drift_gas),
    "Invalid named arguments in GCIMSSample()"
  )
})

test_that("A transposed (wrong orientation) intensity matrix errors immediately", {
  # data should be [drift_time x retention_time]. Giving it transposed
  # ([retention_time x drift_time], as issue #28 did) must fail validation
  # right away, rather than being accepted and misbehaving later on.
  dt <- 1:2
  rt <- 1:3
  wrong_orientation_data <- matrix(seq(length(dt) * length(rt)), nrow = length(rt), ncol = length(dt))
  expect_error(
    GCIMSSample(drift_time = dt, retention_time = rt, data = wrong_orientation_data),
    "data should have as many rows as the drift_time length"
  )
})

test_that("A non-square intensity matrix round-trips correctly through intensity()", {
  dt <- seq(0, 0.4, by = 0.1) # length 5
  rt <- seq(0, 4, by = 0.8)   # length 6
  data <- matrix(seq(length(dt) * length(rt)), nrow = length(dt), ncol = length(rt))
  obj <- GCIMSSample(drift_time = dt, retention_time = rt, data = data)
  expect_equal(dim(intensity(obj)), c(length(dt), length(rt)))
  expect_equal(unname(intensity(obj)), data)
})

