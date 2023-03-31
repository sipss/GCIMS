phenos_to_string <- function(phenotypes) {
  phenotypes <- as.data.frame(phenotypes)
  only_phenotypes <- phenotypes[,c(-1, -2)]  # Remove SampleID and FileName
  num_of_phenotypes <- length(colnames(only_phenotypes))
  if (num_of_phenotypes > 10) {
    to_show <- 5
    phen <- paste0(utils::head(colnames(only_phenotypes), n = to_show), collapse = ", ")
    and_n_more <- paste0(" and ", num_of_phenotypes - to_show, " more phenotypes")
    return(paste0(phen, and_n_more))
  } else if (num_of_phenotypes == 0) {
    return("No phenotypes")
  }
  to_show <- num_of_phenotypes
  phen <- paste0(utils::head(colnames(only_phenotypes), n = to_show), collapse = ", ")
  phen
}
