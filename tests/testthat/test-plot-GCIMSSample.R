gauss <- function(x, center, height, sd) height * exp(-(x - center)^2 / (2 * sd^2))

make_sample <- function() {
  dt <- seq(0, 4, by = 0.02) # 201 pts
  rt <- seq(0, 50, by = 0.5) # 101 pts
  int_mat <- matrix(50, nrow = length(dt), ncol = length(rt)) +
    outer(gauss(dt, 2, 500, 0.1), gauss(rt, 25, 1, 3))
  GCIMSSample(drift_time = dt, retention_time = rt, data = int_mat)
}

test_that("cubic_root_trans() transforms/inverts values symmetrically, including negatives", {
  tr <- cubic_root_trans()

  expect_equal(tr$name, "cubic_root")
  expect_equal(tr$transform(8), 2)
  expect_equal(tr$transform(-8), -2)
  expect_equal(tr$inverse(tr$transform(c(-27, -8, 0, 8, 27))), c(-27, -8, 0, 8, 27))
})

test_that("cubic_root_trans()'s breaks() handles a normal range and all-non-finite input", {
  tr <- cubic_root_trans()

  b <- tr$breaks(c(-8, 8))
  expect_true(length(b) > 0)
  expect_true(all(is.finite(b)))

  expect_equal(tr$breaks(c(NA, NaN, Inf)), numeric(0))
})

test_that("as.data.frame() on a GCIMSSample melts the intensity matrix into long format", {
  s <- GCIMSSample(drift_time = 1:2, retention_time = 1:3, data = matrix(1:6, nrow = 2, ncol = 3))

  df <- as.data.frame(s)

  expect_named(df, c("dt_ms", "rt_s", "Intensity"))
  expect_equal(nrow(df), 6L)
  expect_equal(df$Intensity, 1:6)
})

test_that("as.data.frame() on a GCIMSSample respects dt_range/rt_range and row.names", {
  s <- GCIMSSample(drift_time = 1:2, retention_time = 1:3, data = matrix(1:6, nrow = 2, ncol = 3))

  df <- as.data.frame(s, dt_range = c(1, 1))
  expect_equal(nrow(df), 3L)

  df_rn <- as.data.frame(s, row.names = paste0("r", 1:6))
  expect_equal(rownames(df_rn), paste0("r", 1:6))
})

test_that("plot() on a GCIMSSample returns a ggplot with the expected labels and axis limits", {
  s <- make_sample()

  p <- plot(s)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Drift time (ms)")
  expect_equal(p$labels$y, "Retention time (s)")
  expect_equal(p$labels$fill, "Intensity (a.u.)")
})

test_that("plot() on a GCIMSSample restricts the plotted range to dt_range/rt_range", {
  s <- make_sample()

  p_full <- plot(s)
  p_sub <- plot(s, dt_range = c(1, 2), rt_range = c(10, 20))

  b_full <- ggplot2::ggplot_build(p_full)
  b_sub <- ggplot2::ggplot_build(p_sub)

  full_x_range <- diff(b_full$layout$panel_params[[1]]$x.range)
  sub_x_range <- diff(b_sub$layout$panel_params[[1]]$x.range)
  expect_lt(sub_x_range, full_x_range)
})

test_that("plot(remove_baseline = TRUE) on a GCIMSSample does not error once a baseline is estimated", {
  s <- make_sample()
  s_b <- estimateBaseline(s, dt_peak_fwhm_ms = 0.3, dt_region_multiplier = 6, rt_length_s = 10, remove = FALSE)

  p <- plot(s_b, remove_baseline = TRUE)

  expect_s3_class(p, "ggplot")
})

test_that("plot(intensity_range = ) overrides the auto-computed color scale limits", {
  s <- GCIMSSample(drift_time = 1:5, retention_time = 1:5, data = matrix(1:25, nrow = 5))

  p_auto <- plot(s)
  p_fixed <- plot(s, intensity_range = c(0, 200))

  tr <- cubic_root_trans()
  expect_equal(p_auto$scales$get_scales("fill")$get_limits(), tr$transform(c(1, 25)))
  expect_equal(p_fixed$scales$get_scales("fill")$get_limits(), tr$transform(c(0, 200)))
})

test_that("plot(intensity_range = 'ranged') is the default and scales to the cropped view", {
  s <- GCIMSSample(drift_time = 1:10, retention_time = 1:10, data = matrix(1:100, nrow = 10))

  p_cropped <- plot(s, dt_range = c(1, 1))

  tr <- cubic_root_trans()
  expect_equal(p_cropped$scales$get_scales("fill")$get_limits(), tr$transform(c(1, 91)))
})

test_that("plot(intensity_range = 'global') scales to the sample's full range, ignoring dt_range/rt_range", {
  s <- GCIMSSample(drift_time = 1:10, retention_time = 1:10, data = matrix(1:100, nrow = 10))

  p_global <- plot(s, dt_range = c(1, 1), intensity_range = "global")

  tr <- cubic_root_trans()
  expect_equal(p_global$scales$get_scales("fill")$get_limits(), tr$transform(c(1, 100)))
})

test_that("plot(intensity_range = 'global') errors clearly with remove_baseline = TRUE", {
  s <- GCIMSSample(drift_time = 1:5, retention_time = 1:5, data = matrix(1:25, nrow = 5))
  s <- estimateBaseline(s, dt_peak_fwhm_ms = 1, dt_region_multiplier = 2, rt_length_s = 2)

  expect_error(
    plot(s, remove_baseline = TRUE, intensity_range = "global"),
    "global.*remove_baseline"
  )
})

test_that("plot(intensity_range = list(min=, max=)) resolves each endpoint independently", {
  s <- GCIMSSample(drift_time = 1:10, retention_time = 1:10, data = matrix(1:100, nrow = 10))

  p <- plot(s, dt_range = c(1, 1), intensity_range = list(min = "ranged", max = "global"))

  tr <- cubic_root_trans()
  expect_equal(p$scales$get_scales("fill")$get_limits(), tr$transform(c(1, 100)))
})

test_that("plot() on a GCIMSSample accepts trans as a string, or as a transform object directly", {
  s <- make_sample()

  # Regression test: mat_to_gplot() checked inherits(trans, "trans"), the
  # legacy S3 class name for scales transform objects. scales >= 1.3.0
  # renamed it to "transform", which made any non-string trans (including
  # this package's own cubic_root_trans()) always fail with "unknown trans
  # value".
  expect_s3_class(plot(s, trans = "cubic_root"), "ggplot")
  expect_s3_class(plot(s, trans = "log10"), "ggplot")
  expect_s3_class(plot(s, trans = cubic_root_trans()), "ggplot")
  expect_s3_class(plot(s, trans = scales::identity_trans()), "ggplot")
})

test_that("plot() on a GCIMSSample errors clearly for an invalid trans value", {
  s <- make_sample()

  expect_error(plot(s, trans = "not_a_real_transform"), "unknown trans value")
  expect_error(plot(s, trans = 123), "unknown trans value")
})

test_that("overlay_peaklist() with a single color returns a rect layer, plus a point or blank layer for the apex", {
  peaklist <- data.frame(
    SampleID = c("s1", "s1", "s2"),
    dt_min_ms = c(1, 3, 2), dt_max_ms = c(2, 4, 3),
    rt_min_s = c(1, 3, 2), rt_max_s = c(2, 4, 3)
  )

  out_no_apex <- overlay_peaklist(peaklist)
  expect_length(out_no_apex, 2)
  expect_s3_class(out_no_apex[[1]]$geom, "GeomRect")
  expect_s3_class(out_no_apex[[2]]$geom, "GeomBlank")

  peaklist$dt_apex_ms <- c(1.5, 3.5, 2.5)
  peaklist$rt_apex_s <- c(1.5, 3.5, 2.5)
  out_apex <- overlay_peaklist(peaklist)
  expect_length(out_apex, 2)
  expect_s3_class(out_apex[[1]]$geom, "GeomRect")
  expect_s3_class(out_apex[[2]]$geom, "GeomPoint")
})

test_that("overlay_peaklist() accepts a literal color name for color_by", {
  peaklist <- data.frame(
    SampleID = "s1", dt_min_ms = 1, dt_max_ms = 2, rt_min_s = 1, rt_max_s = 2
  )

  out <- overlay_peaklist(peaklist, color_by = "blue")

  expect_length(out, 2)
})

test_that("overlay_peaklist() with color_by as a column adds a color scale and, past 10 groups, hides the legend", {
  peaklist <- data.frame(
    SampleID = c("s1", "s2"),
    dt_min_ms = c(1, 2), dt_max_ms = c(2, 3), rt_min_s = c(1, 2), rt_max_s = c(2, 3)
  )

  out <- overlay_peaklist(peaklist, color_by = "SampleID")
  geom_classes <- purrr::map_chr(out, function(l) if (inherits(l, "Layer")) class(l$geom)[1] else class(l)[1])
  expect_equal(geom_classes, c("GeomRect", "GeomBlank", "GeomBlank", "ScaleDiscrete"))

  many <- data.frame(
    SampleID = paste0("s", 1:12),
    dt_min_ms = 1:12, dt_max_ms = 2:13, rt_min_s = 1:12, rt_max_s = 2:13
  )
  out_many <- overlay_peaklist(many, color_by = "SampleID")
  geom_classes_many <- purrr::map_chr(out_many, function(l) if (inherits(l, "Layer")) class(l$geom)[1] else class(l)[1])
  expect_equal(geom_classes_many, c("GeomRect", "GeomBlank", "GeomBlank", "ScaleDiscrete", "Guides"))
})

test_that("overlay_peaklist() with color_by as a column AND apex columns adds a point layer too", {
  peaklist <- data.frame(
    SampleID = c("s1", "s2"),
    dt_min_ms = c(1, 2), dt_max_ms = c(2, 3), rt_min_s = c(1, 2), rt_max_s = c(2, 3),
    dt_apex_ms = c(1.5, 2.5), rt_apex_s = c(1.5, 2.5)
  )

  out <- overlay_peaklist(peaklist, color_by = "SampleID")
  geom_classes <- purrr::map_chr(out, function(l) if (inherits(l, "Layer")) class(l$geom)[1] else class(l)[1])

  expect_equal(geom_classes, c("GeomRect", "GeomBlank", "GeomPoint", "ScaleDiscrete"))
})

test_that("overlay_peaklist() recycles a palette smaller than the number of groups", {
  many <- data.frame(
    SampleID = paste0("s", 1:5),
    dt_min_ms = 1:5, dt_max_ms = 2:6, rt_min_s = 1:5, rt_max_s = 2:6
  )

  expect_no_error(
    overlay_peaklist(many, color_by = "SampleID", palette = c("red", "blue"))
  )
})

test_that("overlay_peaklist() coerces a non-data.frame peaklist (e.g. S4Vectors::DataFrame)", {
  peaklist_df <- S4Vectors::DataFrame(SampleID = "s1", dt_min_ms = 1, dt_max_ms = 2, rt_min_s = 1, rt_max_s = 2)

  out <- overlay_peaklist(peaklist_df)

  expect_length(out, 2)
})

test_that("overlay_peaklist() merges in pdata by SampleID when given", {
  peaklist <- data.frame(
    SampleID = c("s1", "s2"),
    dt_min_ms = c(1, 2), dt_max_ms = c(2, 3), rt_min_s = c(1, 2), rt_max_s = c(2, 3)
  )
  pdata <- data.frame(SampleID = c("s1", "s2"), Group = c("A", "B"))

  out <- overlay_peaklist(peaklist, pdata = pdata, color_by = "Group")

  expect_length(out, 4)
})

test_that("mat_to_gplot() derives dt/rt limits from the matrix dimnames when not given explicitly", {
  m <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("1", "2"), c("10", "20", "30")))

  p <- mat_to_gplot(m)

  b <- ggplot2::ggplot_build(p)
  expect_equal(b$layout$panel_params[[1]]$x.range, c(0.95, 2.05))
  expect_equal(b$layout$panel_params[[1]]$y.range, c(9, 31))
})

test_that("overlay_peaklist() validates color_by and mapping_roi", {
  peaklist <- data.frame(
    SampleID = "s1", dt_min_ms = 1, dt_max_ms = 2, rt_min_s = 1, rt_max_s = 2
  )

  expect_error(overlay_peaklist(peaklist, color_by = "not_a_col_or_color"), "neither a valid color name")
  expect_error(overlay_peaklist(peaklist, color_by = 5), "should be .* a string")
  expect_error(
    overlay_peaklist(peaklist, mapping_roi = c(a = "a", b = "b", c = "c", d = "d")),
    "mapping_roi should be a named vector"
  )
})
