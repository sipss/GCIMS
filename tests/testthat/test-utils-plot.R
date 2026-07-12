test_that("mat_to_nativeRaster() encodes a matrix into a nativeRaster with the right shape", {
  colormap <- farver::encode_native(c("#000000", "#333333", "#666666", "#999999", "#CCCCCC", "#FFFFFF"))
  m <- matrix(1:6, nrow = 2, ncol = 3)

  nr <- mat_to_nativeRaster(m, colormap)

  expect_s3_class(nr, "nativeRaster")
  expect_equal(dim(nr), c(ncol(m), nrow(m)))
  expect_equal(attr(nr, "channels"), 4L)
  # more than one color of the map should show up for a varied matrix:
  expect_gt(length(unique(as.vector(unclass(nr)))), 1L)
  expect_true(all(unique(as.vector(unclass(nr))) %in% colormap))
})

test_that("mat_to_nativeRaster() is deterministic", {
  colormap <- farver::encode_native(c("#000000", "#FFFFFF"))
  m <- matrix(1:6, nrow = 2, ncol = 3)

  expect_identical(mat_to_nativeRaster(m, colormap), mat_to_nativeRaster(m, colormap))
})

test_that("mat_to_nativeRaster() accepts a character colormap, giving the same result as pre-encoded colors", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  chr_colors <- c("#000000", "#333333", "#666666", "#999999", "#CCCCCC", "#FFFFFF")

  nr_chr <- mat_to_nativeRaster(m, chr_colors)
  nr_num <- mat_to_nativeRaster(m, farver::encode_native(chr_colors))

  expect_identical(nr_chr, nr_num)
})

test_that("mat_to_nativeRaster() handles a zero-length matrix", {
  colormap <- farver::encode_native(c("#000000", "#FFFFFF"))

  nr <- mat_to_nativeRaster(matrix(numeric(0), nrow = 0, ncol = 0), colormap)

  expect_s3_class(nr, "nativeRaster")
  expect_equal(dim(nr), c(0L, 0L))
  expect_equal(attr(nr, "channels"), 4L)
})

test_that("plot_interactive() wraps a ggplot into a plotly widget", {
  d <- data.frame(x = c(1, 2), y = c(1, 2))
  plt <- ggplot2::ggplot(d) + ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y))

  out <- plot_interactive(plt)

  expect_s3_class(out, "plotly")
})
