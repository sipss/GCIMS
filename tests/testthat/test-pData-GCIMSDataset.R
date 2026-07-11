make_dataset <- function() {
  sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
  pdat <- data.frame(SampleID = "s1", FileName = basename(sample_file))
  GCIMSDataset$new(pData = pdat, base_dir = dirname(sample_file), on_ram = TRUE)
}

test_that("pData() returns the phenotype data as a tibble", {
  ds <- make_dataset()

  out <- pData(ds)

  expect_s3_class(out, "tbl_df")
  expect_equal(as.character(out$SampleID), "s1")
})

test_that("pData<- replaces the phenotype data, keeping extra annotation columns", {
  ds <- make_dataset()
  newp <- pData(ds)
  newp$Group <- "A"

  pData(ds) <- newp

  expect_equal(pData(ds)$Group, "A")
})

test_that("pData<- errors when the number of rows does not match", {
  ds <- make_dataset()

  expect_error(
    {
      pData(ds) <- data.frame(SampleID = c("s1", "s2"), FileName = c("a", "b"))
    },
    "Number of rows expected"
  )
})

test_that("pData<- warns when FileName is changed", {
  ds <- make_dataset()
  newp <- pData(ds)
  newp$FileName <- "different.mea.gz"

  expect_warning(
    {
      pData(ds) <- newp
    },
    "FileName changed"
  )
})
