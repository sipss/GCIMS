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
  require_pkgs(c("httr", "curl"))
  req <- httr::GET("https://api.github.com/repos/sipss/datasets-gcims-demo/contents/2021-mixture-six-ketones-demo")
  httr::stop_for_status(req)
  req_content <- httr::content(req)

  urls <- purrr::map_chr(req_content, "download_url")
  names(urls) <- purrr::map_chr(req_content, "name")

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  purrr::imap(urls, function(url, name) {
    outfile <- file.path(outdir, name)
    if (file.exists(outfile)) {
      return()
    }
    curl::curl_download(url, outfile)
  })
}

