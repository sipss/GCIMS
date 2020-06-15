
readGCIMS <- function(x) {
  library(readr)
  library(sjmisc)
  if (str_contains(x, ".csv")){ ## I ask if my file (x) contains ".csv in the name. It it is a csv, I read the file.
    file <- read_delim(x, ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
    metadata <- file[1:3, ]
    data <- file[-(1:3), ]
    colnames(data) <- data[1, ]
    data <- data[-1, ]
  } else { ## If my file i is not a csv, I show the warning
   warning("This fild is not a .csv file")
  }
  data <- list(metadata, data)
}


