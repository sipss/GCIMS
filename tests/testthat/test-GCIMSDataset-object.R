create_dummy_dataset <- function() {
  s1 <- GCIMSSample(1:3, 1:2, matrix(1, nrow = 3, ncol = 2))
  s2 <- GCIMSSample(1:3, 1:2, matrix(2, nrow = 3, ncol = 2))
  d1 <- GCIMSDataset$new_from_list(
    samples = list(s1 = s1, s2 = s2),
    on_ram = TRUE
  )
  list(s1 = s1, s2 = s2, d1 = d1)
}

test_that("Create GCIMSDataset from list", {
  a <- create_dummy_dataset()
  expect_equal(a$d1$sampleNames, c("s1", "s2"))
})


test_that("Create and copy a dataset", {
  a <- create_dummy_dataset()
  a$d1$userData$mydata <- "potato"
  d2 <- a$d1$copy()
  expect_equal(a$d1$userData$mydata, "potato")
  expect_equal(d2$userData$mydata, "potato")
  d2$userData$mydata <- NULL
  expect_equal(a$d1$userData$mydata, "potato")
})

test_that("Get a sample", {
  a <- create_dummy_dataset()
  s1 <- a$d1$getSample(1)
  s1b <- a$d1$getSample("s1")
  description(a$s1) <- "s1"
  expect_equal(s1, a$s1)
  expect_equal(s1, s1b)
})

