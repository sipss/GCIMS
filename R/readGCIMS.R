#' GC-IMS Samples Reading


#' @param x               Name of the sample
#' @return An R object that contains the matrix and the retention times
#' @family Reading functions
#' @export
#' @importFrom readr read_csv
#' @importFrom sjmisc str_contains
#' @importFrom R.matlab readMat
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

readGCIMS <- function(dir_in, dir_out) {
  setwd(dir_in)
  metadata <- NULL
  files <- list.files(pattern = ".csv")
  for (i in 1:length(files)){
    if (str_contains(files[i], ".csv")){ ## I ask if my file (x) contains ".csv in the name. It it is a csv, I read the file.
      file <- readr::read_csv(files[i], skip = 130)
      dd <- file[-1, -c(1:2)]
      dd <- list(metadata, dd)
      setwd(dir_out)
      saveRDS(data, file = paste0("M", i, ".rds"))
      setwd(dir_in)
    }else if (sjmisc::str_contains(x, ".mat")) {
      dd <- R.matlab::readMat(x)[[1]]
    }else { ## If my file i is not a csv, I show the warning
      warning("This fild is not a .csv or a .mat file")
    }
  }
}
