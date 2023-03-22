test_that("plotRIC returns a gg object", {
  s1 <- GCIMSSample(1:3, 1:5, matrix(1:15, nrow = 3))
  ds <- GCIMSDataset_fromList(list(s1 = s1, s2 = s1), on_ram = FALSE)
  gplt <- plotRIC(ds)
  expect_s3_class(gplt, "gg")
})
