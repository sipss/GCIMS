#' Align a [GCIMSDataset] using the POW algorithm
#'
#' Uses the POW algorithm to generate waprs and align the samples
#'
#' To use this function you need to install the "pow" package hosted on github, to do this run he next line:
#' remotes::install_github("sipss/pow"),
#' also you will need the "pracma" package hosted in CRAN
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @param lambdas vector of lambdas to be tested for the POW algorithm
#' @param p points to be taken for validation for the POW algorithm
#' @param max_it maximum iterations for applying the POW algorithm
#' @return The modified [GCIMSDataset]
#' @export

# methods::setMethod(
#   "align_pow",
#   "GCIMSDataset",
#   function(object,
#            lambdas = pracma::logspace(-2, 4, 31),
#            p = 10,
#            max_it = 1000){
#     rics <- object$getRIC()
#     idx_y <- pow::select_reference(rics)
#     X <- as.matrix(rics[-idx_y,])
#     y <- as.matrix(rics[idx_y,])
#     n_samples <- nrow(X)
#     m <- length(y)
#     v <- rep(1,m)
#     iv <- seq(2, m - 1, by = p)
#     v[iv] <- 0
#     W <- Matrix::Diagonal(x = v)
#     val<-pow::compute_val_error(X, y, W, iv, lambdas)
#     opar <- pow::optimize_params(n_samples, lambdas, val$ti_ix, val$e_ix)
#     best_lambdas2 <- opar$best_params
#
#     w<-mapply(pow::pow,x=as.list(data.frame(t(X))), lambda2=as.list(best_lambdas2),MoreArgs = list(y=y,max_it = max_it))
#     w<-as.list((as.data.frame(w)))
#     w<-append(w,list(1:m),idx_y-1)
#     names(w)<-object$sampleNames
#
#     tis_matrix <- getTIS(object)
#     dt <- dtime(object)
#     rip_position <- apply(tis_matrix, 1L, which.max)
#     rip_ref_idx <- round(stats::median(rip_position, na.rm = TRUE))
#     rip_ref_ms <- dt[rip_ref_idx]
#
#     delayed_op <- DelayedOperation(
#       name = "align_pow",
#       fun = align_pow,
#       params = list(rip_ref_ms = rip_ref_ms),
#       params_iter = list(w=w)
#     )
#     object$appendDelayedOp(delayed_op)
#     object <- GCIMS:::extract_RIC_and_TIS(object)
#     invisible(object)
#   }
# )
