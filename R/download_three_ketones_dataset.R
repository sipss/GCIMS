#' Download three samples (6-ketone mixture)
#'
#' This function downloads three samples in .mea.gz format. It is useful
#' to run the introductory vignette.
#'
#' @param outdir Name of the directory where the samples will be saved
#'
#' @return Nothing (Files are created in the given folder)
#'
#'
#' @examples
#' \dontrun{
#' download_three_ketones_dataset(outdir = "sample_dataset")
#' list.files("sample_dataset")
#' }
#'
#' @export
download_three_ketones_dataset <- function(outdir = "2021-mixture-six-ketones-demo") {
  require_pkgs("curl")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  mea_files <- c("220221_102919.mea", "220228_114404.mea", "220505_110658.mea")
  full_files <- file.path(outdir, mea_files)
  if (all(file.exists(full_files))) {
    return(full_files)
  }
  url <- "https://zenodo.org/record/7941230/files/AnalyticalStandardsDemo.zip?download=1"
  tmp_zipfile <- tempfile(fileext=".zip")
  curl::curl_download(url, tmp_zipfile)
  utils::unzip(
      tmp_zipfile,
      junkpaths = TRUE, exdir = outdir
  )
  full_files
}

