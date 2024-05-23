#' Create a table for defining dataset annotations
#'
#' To process an entire dataset, we need a table that describes the samples,
#' and you may want to add for the analysis.
#'
#' The table needs to have at least two columns, one with the file name of
#' the sample (`FileName`) and another one with the sample name (`SampleID`), that you can set
#' as you like. Besides, you can add additional columns with any
#' metadata/annotations/phenotypes you may consider relevant.
#'
#' This function will help you list all the samples from a directory to a table. The example
#' below will show you how to save this table as an Excel or CSV file, for you
#' to conveniently modify it and how you can read it back for further analysis.
#'
#' @param samples_dir A directory that contains samples
#' @param glob One or more globs for sample extensions. See the examples.
#' @param recursive Whether to look for samples into `samples_dir` subdirectories
#' @param verbose If set to `TRUE` it prints instructions
#'
#' @return A data frame with the SampleID and FileName columns
#' @export
#'
#' @examples
#' # How to create the annotations table:
#' #
#' # First you must tell R where your samples are. Please change "samples_dir"
#' # below to your samples directory. On Windows you can use:
#' #     samples_dir <- choose.dir(getwd(), "Choose the folder where the samples are")
#' # On other systems you can use:
#' #     library(tcltk)
#' #     samples_dir <- tclvalue(tkchooseDirectory())
#' # In this example we use a folder with some demo files:
#' samples_dir <- system.file("extdata", "sample_formats", package = "GCIMS")
#' # Then you need to provide the extension to look into. If you use `glob = "*.*"` you
#' # will catch all files and you can filter the annotations table afterwards:
#' annotations <- create_annotations_table(samples_dir, glob = "*.mea.gz")
#' # You can write the annotations table to an Excel or a CSV file:
#' # For Excel you may need to install the writexl package:
#' #    install.packages("writexl")
#' # And then you can use:
#' #    writexl::write_xlsx(annotations, "annotations.xlsx")
#' # For csv just use:
#' #    write.csv(annotations, "annotations.csv")
#' #
#' # Modify manually the excel or CSV file
#' #
#' # Read it again into R as follows:
#' #
#' # For Excel you may need to install the readxl package:
#' #    install.packages("readxl")
#' # And then you can use:
#' #    annotations <- readxl::read_excel("annotations.xlsx")
#' # For csv just use:
#' #    annotations <- read.csv("annotations.csv")
create_annotations_table <- function(
    samples_dir,
    glob = c("*.mea", "*.mea.gz"),
    recursive = TRUE,
    verbose = TRUE
) {
  filenames <- list.files(
    samples_dir,
    pattern = paste0(utils::glob2rx(glob), collapse = "|"),
    recursive = recursive,
    include.dirs = FALSE
  )
  annotations <- tibble::tibble(
    SampleID = sapply(filenames, function(x) sub("\\.mea(\\.gz)?$", "", basename(x))),
    FileName = filenames
  )
  if (verbose) {
    cli_inform(
      message = c(
        "An annotation table was created",
        "i" = "The table now includes {nrow(annotations)} samples",
        "i" = "Feel free to edit the table to include additional annotations as extra columns (for groups, phenotypes...), please avoid renaming the SampleID and FileName columns",
        "i" = "You may freely add/remove rows to include/exclude additional samples",
        "i" = paste0(
          "For editing, you can do it from R or you can check the examples at {.code help(\"create_annotations_table\")} to",
          "learn how to save the table to an Excel file at your convenience"
        )
      )
    )
  }
  annotations
}
