test_that("sample_name_or_number_to_both() resolves a single numeric index", {
  out <- sample_name_or_number_to_both(1, c("a", "b", "c"))

  expect_equal(out$idx, 1)
  expect_equal(out$name, "a")
  expect_equal(out$logical, c(TRUE, FALSE, FALSE))
})

test_that("sample_name_or_number_to_both() errors when a numeric index is out of range", {
  expect_error(
    sample_name_or_number_to_both(5, c("a", "b", "c")),
    "between 1 and 3"
  )
})

test_that("sample_name_or_number_to_both() resolves sample names", {
  out <- sample_name_or_number_to_both(c("a", "c"), c("a", "b", "c"))

  expect_equal(out$idx, c(1, 3))
  expect_equal(out$logical, c(TRUE, FALSE, TRUE))
})

test_that("sample_name_or_number_to_both() errors listing sample names that don't exist", {
  expect_error(
    sample_name_or_number_to_both(c("a", "z"), c("a", "b", "c")),
    "Missing sample names: z"
  )
})

test_that("sample_name_or_number_to_both() resolves a logical vector", {
  out <- sample_name_or_number_to_both(c(TRUE, FALSE, TRUE), c("a", "b", "c"))

  expect_equal(out$idx, c(1, 3))
  expect_equal(out$name, c("a", "c"))
  expect_equal(out$logical, c(TRUE, FALSE, TRUE))
})

test_that("sample_name_or_number_to_both() treats NA as FALSE in a logical vector", {
  out <- sample_name_or_number_to_both(c(TRUE, NA, TRUE), c("a", "b", "c"))

  expect_equal(out$idx, c(1, 3))
  expect_equal(out$logical, c(TRUE, FALSE, TRUE))
})

test_that("sample_name_or_number_to_both() errors when the logical vector has the wrong length", {
  expect_error(
    sample_name_or_number_to_both(c(TRUE, FALSE), c("a", "b", "c")),
    "length 3"
  )
})

test_that("sample_name_or_number_to_both() errors for an unsupported type", {
  expect_error(
    sample_name_or_number_to_both(list(1), c("a", "b", "c")),
    "numeric vector, a character vector or a logical vector"
  )
})
