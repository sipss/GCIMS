quantile_criteria <- function(noise_quantile = 0.25) {
  force(noise_quantile)
  setup <- function(intens, half_min_size) {
    thr <- stats::quantile(intens, noise_quantile)
    list(
      binarized = intens <= thr
    )
  }
  eval_region <- function(intens, region, params) {
    all(!params$binarized[
      region["row_min"]:region["row_max"],
      region["col_min"]:region["col_max"]
    ])
  }
  list(
    name = "quantile_criteria",
    setup = setup,
    eval_region = eval_region
  )
}

no_positive_outliers <- function(k = 3) {
  setup <- function(intens, half_min_size) {
    list()
  }
  eval_region <- function(intens, region, params) {
    intens_region <- intens[
      region["row_min"]:region["row_max"],
      region["col_min"]:region["col_max"]
    ]
    threshold <- stats::median(intens_region) + k*stats::mad(intens_region)
    all(intens_region < threshold)
  }
  list(
    name = "no_positive_outliers",
    setup = setup,
    eval_region = eval_region
  )
}

# no_positive_outliers <- function(k = 3) {
#   force(k)
#   setup <- function(intens, half_min_size) {
#     thr <- stats::median(intens) + k*stats::mad(intens)
#     list(
#       binarized = intens <= thr
#     )
#   }
#   eval_region <- function(intens, region, params) {
#     all(!params$binarized[
#       region["row_min"]:region["row_max"],
#       region["col_min"]:region["col_max"]
#     ])
#   }
#   list(
#     name = "no_positive_outliers",
#     setup = setup,
#     eval_region = eval_region
#   )
# }

find_all_regions_without_peaks <- function(intens, half_min_size = c(20, 20), criteria = quantile_criteria(noise_quantile = 0.25), find_all = FALSE) {
  params <- criteria$setup(intens, half_min_size)

  center_i <- seq.int(from = 1L+half_min_size[1], to = nrow(intens) - half_min_size[1], by = 2L*half_min_size[1]+1L)
  center_j <- seq.int(from = 1L+half_min_size[2], to = ncol(intens) - half_min_size[2], by = 2L*half_min_size[2]+1L)

  out <- list()
  for (i in center_i) {
    for (j in center_j) {
      region <- c(
        "row_min" = i - half_min_size[1],
        "row_max" = i + half_min_size[1],
        "col_min" = j - half_min_size[2],
        "col_max" = j + half_min_size[2]
      )
      if (criteria$eval_region(intens, region, params)) {
        if (!find_all) {
          return(region)
        }
        out <- c(out, list(region))
      }
    }
  }
  if (!find_all) {
    stop("Could not find a region without peaks")
  }
  df <- dplyr::bind_rows(!!!out)
  df$criteria <- criteria$name
  df
}

infile <- "/home/sergio/sync_ibec/Escritorio/work/ibec/src/GCIMS-Samples/2022-for-peak-clustering/alignr/M1.rds"
sample <- readRDS(infile)

intens <- sample$data$data_df


rectangulos <- dplyr::bind_rows(

  find_all_regions_without_peaks(
    intens,
    half_min_size = c(10, 10),
    criteria = quantile_criteria(noise_quantile = 0.25),
    find_all = TRUE
  ),

  find_all_regions_without_peaks(
    intens,
    half_min_size = c(10, 10),
    criteria = no_positive_outliers(k=3),
    find_all = TRUE
  )
)

intens_long <- data.frame(
  drift_time = rep(sample$data$drift_time, times = length(sample$data$retention_time)),
  retention_time = rep(sample$data$retention_time, each = length(sample$data$drift_time)),
  intensity = as.numeric(intens)
)

intens_long$pt <- seq_len(nrow(intens_long))



library(ggplot2)

plt1 <- ggplot() +
  geom_contour_filled(
    mapping = aes(x = drift_time, y = retention_time, z = intensity^(1/3)),
    data = dplyr::filter(intens_long, pt %% 100 == 0),
    show.legend = FALSE
  )
plt2 <- ggplot() +
  geom_rect(
    aes(
      xmin = sample$data$drift_time[row_min],
      xmax = sample$data$drift_time[row_max],
      ymin = sample$data$retention_time[col_min],
      ymax = sample$data$retention_time[col_max],
      fill = criteria
    ),
    data = rectangulos
  ) +
  facet_wrap(~criteria)

cowplot::plot_grid(plt1,
                   plt2 + guides(fill="none"),
                   rel_widths = c(1, 2))
