test_that("DelayedOperation validates its name", {
  op <- DelayedOperation(name = "my-op1")
  expect_equal(name(op), "my-op1")

  expect_error(
    DelayedOperation(name = "1 bad name!"),
    "not a valid operation name"
  )
})

test_that("modifiesSample is TRUE only when fun is provided", {
  op_with_fun <- DelayedOperation(name = "withfun", fun = function(x) x + 1)
  op_without_fun <- DelayedOperation(name = "nofun")

  expect_true(modifiesSample(op_with_fun))
  expect_false(modifiesSample(op_without_fun))
})

test_that("apply_op_to_sample runs fun with params and reports needs_resaving", {
  op <- DelayedOperation(name = "add", fun = function(x, n) x + n, params = list(n = 5))

  out <- apply_op_to_sample(op, sample = 10, sample_name = "s1")

  expect_equal(out$sample, 15)
  expect_true(out$needs_resaving)
  expect_null(out$extracted_obj)
})

test_that("apply_op_to_sample with only fun_extract leaves the sample untouched", {
  op <- DelayedOperation(name = "extract", fun_extract = function(x) x * 2)

  out <- apply_op_to_sample(op, sample = 10, sample_name = "s1")

  expect_equal(out$sample, 10)
  expect_false(out$needs_resaving)
  expect_equal(out$extracted_obj, 20)
})

test_that("apply_op_to_sample extracts from the sample AFTER fun modifies it", {
  op <- DelayedOperation(
    name = "both",
    fun = function(x, n) x + n, params = list(n = 5),
    fun_extract = function(x) x * 10
  )

  out <- apply_op_to_sample(op, sample = 10, sample_name = "s1")

  expect_equal(out$sample, 15)
  expect_equal(out$extracted_obj, 150) # 150 = 15 * 10, not 10 * 10
})

test_that("apply_op_to_sample passes params_iter values keyed by sample name", {
  op <- DelayedOperation(
    name = "iter",
    fun = function(x, n) x + n,
    params_iter = list(n = list(s1 = 1, s2 = 100))
  )

  out_s1 <- apply_op_to_sample(op, sample = 10, sample_name = "s1")
  out_s2 <- apply_op_to_sample(op, sample = 10, sample_name = "s2")

  expect_equal(out_s1$sample, 11)
  expect_equal(out_s2$sample, 110)
})

test_that("describeAsList returns just the name when there are no params", {
  op <- DelayedOperation(name = "simple-op")

  expect_equal(describeAsList(op), "simple-op")
})

test_that("describeAsList nests params under the operation name", {
  op <- DelayedOperation(name = "with-params", params = list(a = 1, b = "hello"))

  out <- describeAsList(op)

  expect_named(out, "with-params")
  expect_equal(out[["with-params"]], list(a = 1, b = "hello"))
})

test_that("describeAsList truncates long atomic vectors and summarizes data.frame/function params", {
  op <- DelayedOperation(
    name = "big-params",
    params = list(
      v = 1:20,
      df = data.frame(x = 1:3, y = 4:6),
      f = mean
    )
  )

  out <- describeAsList(op)[["big-params"]]

  expect_match(out$v, "numeric vector of 20 elements")
  expect_match(out$df, "data.frame with 3 rows and 2 columns")
  expect_equal(out$f, "< function >")
})

test_that("hashableDelayedOp is reproducible across structurally-identical closures", {
  op_a <- DelayedOperation(name = "h", fun = function(x) x + 1, params = list(n = 1))
  op_b <- DelayedOperation(name = "h", fun = function(x) x + 1, params = list(n = 1))

  expect_identical(hashableDelayedOp(op_a), hashableDelayedOp(op_b))
})

test_that("show.DelayedOperation prints without erroring", {
  op <- DelayedOperation(name = "printable", params = list(a = 1))

  expect_output(print(op), "printable")
})
