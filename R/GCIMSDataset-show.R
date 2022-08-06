phenos_to_string <- function(object) {
  phenotypes <- pData(object)
  only_phenotypes <- phenotypes[,c(-1, -2)]
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

methods::setMethod(
  "describeAsList", "GCIMSDataset",
  function(object) {
    out <- list()
    root_txt <- "A GCIMSDataset with"
    out[[root_txt]] <- list()
    out[[root_txt]][[1]] <- paste0(length(sampleNames(object)), " samples")
    # Phenotype:
    out[[root_txt]][[2]] <- phenos_to_string(object)
    # Previous operations
    pops <- object@envir$previous_ops
    pops <- purrr::keep(pops, modifiesSample)
    pops <- purrr::map(pops, describeAsList)
    if (length(pops) > 0) {
      out[[root_txt]][[3]] <- list("History" = pops)
    } else {
      out[[root_txt]][[3]] <- "No previous history"
    }
    # Pending operations
    pops <- object@envir$delayed_ops
    pops <- purrr::keep(pops, modifiesSample)
    pops <- purrr::map(pops, describeAsList)
    if (length(pops) > 0) {
      out[[root_txt]][[4]] <- list("Pending operations" = pops)
    } else {
      out[[root_txt]][[4]] <- "No pending operations"
    }
    out
  }
)

setMethod(
  "show",
  "GCIMSDataset",
  function(object) {
    outstring <- yaml::as.yaml(describeAsList(object))
    cat(outstring)
  }
)
