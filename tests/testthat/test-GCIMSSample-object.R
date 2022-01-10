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
    "Invalid named arguments in GCIMSSample initialization"
  )
})

