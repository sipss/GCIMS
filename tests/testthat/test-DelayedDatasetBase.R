test_that("DelayedDatasetBase's virtual methods abort with a clear message", {
  base <- DelayedDatasetBase$new()

  expect_error(base$getSample(1, dataset = NULL), "virtual function")
  expect_error(base$sampleNames, "virtual active binding")

  # realize() short-circuits when there are no queued ops, so queue one to
  # actually reach the virtual realize_impl():
  base$appendDelayedOp(DelayedOperation(name = "op", fun = function(x) x))
  expect_error(base$realize(dataset = NULL), "virtual function")
})

test_that("DelayedDatasetBase tracks queued operations", {
  base <- DelayedDatasetBase$new()
  expect_false(base$hasDelayedOps())

  base$appendDelayedOp(DelayedOperation(name = "op1", fun = function(x) x + 1))

  expect_true(base$hasDelayedOps())
  pending <- base$pending_as_list()
  expect_named(pending, "Queued operations")
  expect_equal(base$history_as_list(), "No previous history")
})

test_that("DelayedDatasetRAM realizes queued operations and aggregates extracted results", {
  ds <- DelayedDatasetRAM$new(samples = list(a = 1, b = 2))
  extracted <- NULL
  op <- DelayedOperation(
    name = "add10",
    fun = function(x) x + 10,
    fun_extract = function(x) x,
    fun_aggregate = function(dataset, extracted_result) {
      extracted <<- extracted_result
      dataset
    }
  )
  ds$appendDelayedOp(op)

  ds$realize(dataset = list())

  expect_equal(ds$getSample("a", dataset = list()), 11)
  expect_equal(ds$getSample("b", dataset = list()), 12)
  expect_equal(extracted, list(a = 11, b = 12))
  expect_false(ds$hasDelayedOps())
})

test_that("DelayedDatasetBase registerOptimization is applied before realizing", {
  ds <- DelayedDatasetRAM$new(samples = list(a = 1))
  ds$appendDelayedOp(DelayedOperation(name = "add1", fun = function(x) x + 1))
  ds$appendDelayedOp(DelayedOperation(name = "add100", fun = function(x) x + 100))
  # Drop the first queued operation before it runs:
  ds$registerOptimization(function(ops) ops[-1])

  ds$realize(dataset = list())

  expect_equal(ds$getSample("a", dataset = list()), 101)
})

test_that("DelayedDatasetRAM requires unique, non-empty sample names", {
  expect_error(
    DelayedDatasetRAM$new(samples = list(a = 1, a = 2)),
    "unique and non-empty"
  )
})

test_that("DelayedDatasetRAM$subset keeps only the requested samples", {
  ds <- DelayedDatasetRAM$new(samples = list(a = 1, b = 2, c = 3))

  ds$subset(c("a", "c"))

  expect_equal(ds$sampleNames, c("a", "c"))
})

test_that("DelayedDatasetRAM$subset refuses to run with pending operations", {
  ds <- DelayedDatasetRAM$new(samples = list(a = 1))
  ds$appendDelayedOp(DelayedOperation(name = "op", fun = function(x) x))

  expect_error(ds$subset("a"), "pending operations")
})

test_that("DelayedDatasetRAM$sampleNames setter renames samples and rejects duplicates", {
  ds <- DelayedDatasetRAM$new(samples = list(a = 1, b = 2))

  ds$sampleNames <- c("x", "y")
  expect_equal(ds$sampleNames, c("x", "y"))

  expect_error(
    {
      ds$sampleNames <- c("x", "x")
    },
    "unique and not missing"
  )
})
