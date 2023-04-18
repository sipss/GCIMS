create_dummy_dataset <- function(on_ram = TRUE) {
  s1 <- GCIMSSample(1:3, 1:2, matrix(1, nrow = 3, ncol = 2))
  s2 <- GCIMSSample(1:3, 1:2, matrix(2, nrow = 3, ncol = 2))
  d1 <- GCIMSDataset$new_from_list(
    samples = list(s1 = s1, s2 = s2),
    on_ram = on_ram,
    scratch_dir = NULL
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


test_that("subset on RAM works", {
  a <- create_dummy_dataset(on_ram = TRUE)
  b <- a$d1$subset(samples = "s2", inplace = FALSE)
  expect_equal(a$d1$sampleNames, c("s1", "s2"))
  expect_equal(b$sampleNames, "s2")
})


test_that("subset on disk works", {
  d1 <- create_dummy_dataset(on_ram = FALSE)$d1
  d2 <- d1$subset(samples = "s2", inplace = FALSE, new_scratch_dir = NULL)
  expect_equal(d1$sampleNames, c("s1", "s2"))
  expect_equal(d2$sampleNames, "s2")
  # Different scratch directory:
  expect_true(dirname(d1$getCurrentDir()) != dirname(d2$getCurrentDir()))
  # Same step:
  expect_true(basename(d1$getCurrentDir()) == basename(d2$getCurrentDir()))
  expect_identical(list.files(d1$getCurrentDir()), c("dataset.rds", "sample_s1.rds", "sample_s2.rds"))
  expect_identical(list.files(d2$getCurrentDir()), c("dataset.rds", "sample_s2.rds"))
})
