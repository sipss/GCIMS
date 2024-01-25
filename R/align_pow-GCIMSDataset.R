#' Align a [GCIMSDataset] using the POW algorithm
#'
#' Uses the POW algorithm to generate waprs and align the samples
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @param lambdas vector of lambdas to be tested for the POW algorithm
#' @param p points to be taken for validation for the POW algorithm
#' @param max_it maximum iterations for applying the POW algorithm
#' @return The modified [GCIMSDataset]
#' @export

methods::setMethod(
  "align_pow",
  "GCIMSDataset",
  function(object,
           lambdas = pracma::logspace(-2, 4, 31),
           p = 10,
           max_it = 1000){
    rics <- object$getRIC()
    idx_y <- Align::select_reference(rics)
    X <- as.matrix(rics[-idx_y,])
    y <- as.matrix(rics[idx_y,])
    n_samples <- nrow(X)
    m <- length(y)
    v <- rep(1,m)
    iv <- seq(2, m - 1, by = p)
    v[iv] <- 0
    W <- Matrix::Diagonal(x = v)
    val<-Align::compute_val_error(X,y,W,iv,lambdas)
    opar <- Align::optimize_params(n_samples, lambdas, val$ti_ix, val$e_ix)
    best_lambdas2 <- opar$best_params

    w<-mapply(Align::pow,x=as.list(data.frame(t(X))), lambda2=as.list(best_lambdas2),MoreArgs = list(y=y,max_it = max_it))
    w<-as.list((as.data.frame(w)))
    w<-append(w,list(1:m),idx_y-1)
    names(w)<-object$sampleNames

    delayed_op <- DelayedOperation(
      name = "align_pow",
      fun = align_pow,
      params_iter = list(w=w)
    )
    object$appendDelayedOp(delayed_op)
    object <- GCIMS:::extract_RIC_and_TIS(object)
    invisible(object)
  }
)
