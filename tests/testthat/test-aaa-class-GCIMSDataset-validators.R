make_sample <- function() GCIMSSample(1:2, 1:2, matrix(1, nrow = 2, ncol = 2))

test_that("validate_on_ram accepts TRUE/FALSE-like values and rejects the rest", {
  expect_true(validate_on_ram(TRUE))
  expect_false(validate_on_ram(FALSE))
  expect_true(validate_on_ram(1))

  expect_error(validate_on_ram(c(TRUE, FALSE)), "on_ram must be either TRUE or FALSE")
  expect_error(validate_on_ram(NA), "on_ram must be either TRUE or FALSE")
})

test_that("validate_keep_intermediate requires a single logical", {
  expect_true(validate_keep_intermediate(TRUE))
  expect_error(validate_keep_intermediate("yes"), "must be either TRUE or FALSE")
})

test_that("validate_parser accepts 'default' or a function", {
  expect_equal(validate_parser("default"), "default")
  f <- function(x) x
  expect_identical(validate_parser(f), f)
  expect_error(validate_parser("bogus"), "must be either 'default' or a function")
})

test_that("validate_samples accepts NULL and a named list of GCIMSSample objects", {
  expect_null(validate_samples(NULL))

  s <- make_sample()
  out <- validate_samples(list(a = s))
  expect_named(out, "a")
})

test_that("validate_samples rejects a non-list", {
  expect_error(validate_samples(1:3), "named list of GCIMSSample")
})

test_that("validate_samples auto-names an unnamed list with a warning", {
  s <- make_sample()
  expect_warning(
    out <- validate_samples(list(s)),
    "Using Sample1"
  )
  expect_named(out, "Sample1")
})

test_that("validate_samples rejects duplicated or empty names, and non-GCIMSSample elements", {
  s <- make_sample()
  expect_error(validate_samples(list(a = s, a = s)), "should be unique")
  expect_error(validate_samples(list(a = 1)), "not of type `GCIMSSample`")
})

test_that("validate_pData requires a SampleID column and errors with a clear message otherwise", {
  expect_error(
    validate_pData(data.frame(Foo = 1)),
    "should have a SampleID column"
  )
})

test_that("validate_pData normalizes SampleID/FileName columns and coerces integer IDs", {
  pd <- validate_pData(data.frame(SampleID = c("s1", "s2"), Extra = "x"))
  expect_equal(colnames(pd)[1:2], c("SampleID", "FileName"))
  expect_equal(as.character(pd$SampleID), c("s1", "s2"))

  expect_warning(
    pd2 <- validate_pData(data.frame(SampleID = 1:2)),
    "coerced to a character column"
  )
  expect_equal(as.character(pd2$SampleID), c("1", "2"))
})

test_that("validate_pData_samples defaults to an empty pData/samples pair when both are NULL", {
  out <- validate_pData_samples(NULL, NULL)
  expect_equal(nrow(out$pData), 0L)
  expect_equal(length(out$samples), 0L)
})

test_that("validate_pData_samples builds pData from sample names when pData is NULL", {
  s <- make_sample()
  out <- validate_pData_samples(NULL, list(a = s))
  expect_equal(as.character(out$pData$SampleID), "a")
})

test_that("validate_pData_samples errors with a readable message when samples are missing from pData", {
  s <- make_sample()
  expect_error(
    suppressWarnings(validate_pData_samples(data.frame(SampleID = "b"), list(a = s))),
    "missing in pData"
  )
})

test_that("check_files errors listing missing files", {
  expect_error(
    check_files("does_not_exist.mea", tempdir()),
    "does_not_exist.mea"
  )
})

test_that("check_files passes silently when all files exist", {
  f <- tempfile(fileext = ".mea")
  file.create(f)
  on.exit(unlink(f))

  expect_no_error(check_files(basename(f), dirname(f)))
})

test_that("abort_if_errors aborts only when there are errors", {
  expect_error(abort_if_errors("bad thing", title = "Oops"), "Oops")
  expect_no_error(abort_if_errors(character(0)))
})

test_that("validate_base_dir requires an existing directory and normalizes the path", {
  out <- validate_base_dir(tempdir())
  expect_true(dir.exists(out))
  expect_error(validate_base_dir("/does/not/exist/xyz123"), "does not exist")
})

test_that("validate_scratch_dir creates a temporary directory when NULL", {
  out <- validate_scratch_dir(NULL)
  expect_true(dir.exists(out))
})

test_that("optimize_RIC_TIS moves extract_dtime_rtime/extract_RIC_and_TIS operations to the end", {
  op1 <- DelayedOperation(name = "smooth", fun = function(x) x)
  op2 <- DelayedOperation(name = "extract_dtime_rtime", fun_extract = function(x) x)
  op3 <- DelayedOperation(name = "detectPeaks", fun = function(x) x)
  op4 <- DelayedOperation(name = "extract_RIC_and_TIS", fun_extract = function(x) x)

  out <- optimize_RIC_TIS(list(op1, op2, op3, op4))

  expect_equal(purrr::map_chr(out, name), c("smooth", "detectPeaks", "extract_dtime_rtime", "extract_RIC_and_TIS"))
})

test_that("optimize_RIC_TIS is a no-op when neither extraction operation is queued", {
  op1 <- DelayedOperation(name = "smooth", fun = function(x) x)
  op3 <- DelayedOperation(name = "detectPeaks", fun = function(x) x)

  out <- optimize_RIC_TIS(list(op1, op3))

  expect_identical(purrr::map_chr(out, name), c("smooth", "detectPeaks"))
})

test_that("read_sample dispatches .mea/.mea.gz to read_mea and reads a real sample", {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")

  out <- read_sample(basename(sample_file), base_dir = dirname(sample_file))

  expect_s4_class(out, "GCIMSSample")
})

test_that("read_sample errors on an unsupported extension, naming the actual file", {
  expect_error(
    read_sample("foo.xyz", base_dir = "/some/dir"),
    "/some/dir/foo.xyz"
  )
})

test_that("read_sample dispatches to a custom parser function when given one", {
  out <- read_sample("whatever.xyz", base_dir = tempdir(), parser = function(filename) list(custom = filename))

  expect_named(out, "custom")
})
