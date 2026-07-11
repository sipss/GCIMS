make_samples_dir <- function() {
  d <- tempfile("annot_")
  dir.create(d)
  dir.create(file.path(d, "sub"))
  file.create(file.path(d, "a.mea"))
  file.create(file.path(d, "b.mea.gz"))
  file.create(file.path(d, "c.txt"))
  file.create(file.path(d, "sub", "d.mea"))
  d
}

test_that("create_annotations_table() lists matching files recursively by default, stripping .mea/.mea.gz from SampleID", {
  d <- make_samples_dir()

  out <- create_annotations_table(d, verbose = FALSE)

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("SampleID", "FileName"))
  expect_setequal(out$FileName, c("a.mea", "b.mea.gz", file.path("sub", "d.mea")))
  expect_setequal(out$SampleID, c("a", "b", "d"))
})

test_that("create_annotations_table(recursive = FALSE) only lists files in the top-level directory", {
  d <- make_samples_dir()

  out <- create_annotations_table(d, recursive = FALSE, verbose = FALSE)

  expect_setequal(out$FileName, c("a.mea", "b.mea.gz"))
})

test_that("create_annotations_table() respects a custom glob", {
  d <- make_samples_dir()

  out <- create_annotations_table(d, glob = "*.mea", recursive = FALSE, verbose = FALSE)

  expect_equal(out$FileName, "a.mea")
})

test_that("create_annotations_table(verbose = TRUE) informs how many samples were found", {
  d <- make_samples_dir()

  expect_message(
    out <- create_annotations_table(d, verbose = TRUE),
    "3 samples"
  )
  expect_equal(nrow(out), 3)
})
