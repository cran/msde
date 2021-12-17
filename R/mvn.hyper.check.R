#' Argument checking for the default multivariate normal prior.
#'
#' @param hyper The normal prior's hyperparameters: `NULL`, or a list with elements `mu` and `Sigma`, corresponding to a named mean vector and variance matrix (see Details).
#' @param param.names Vector of parameter names (see Details).
#' @param data.names Vector of data names (see Details).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`mean`}{The mean vector.}
#'   \item{`cholSd`}{The upper upper Cholesky factor of the variance matrix.}
#'   \item{`thetaId`}{The index of the corresponding variables in `theta`.}
#'   \item{`xId`}{The index of the corresponding variables in `x0`.}
#' }
#'
#' @details This function is not meant to be called directly by the user, but rather to parse the hyper-parameters of a default multivariate normal prior distribution to be passed to the C++ code in [sde.prior()] and [sde.post()].  This default prior is multivariate normal on the elements of `(theta, x0)` specified by each of `names(mu)`, `rownames(Sigma)`, and `colnames(Sigma)`.  The remaining components are given Lebesgue priors, or a full Lebesgue prior if `hyper == NULL`.  If the names of `mu` and `Sigma` are inconsistent an error is thrown.
#' @export
mvn.hyper.check <- function(hyper, param.names, data.names) {
  nparams <- length(param.names)
  ndims <- length(data.names)
  var.names <- c(param.names, data.names)
  prior.names <- c("mu", "Sigma")
  if(!is.null(hyper)) {
    # error checking
    if(is.null(names(hyper)) ||
       !identical(sort(names(hyper)), sort(prior.names))) {
      stop("hyper must be NULL or a list with elements mu and Sigma.")
    }
    # check argument names
    mu.names <- names(hyper$mu)
    Sigma.rnames <- rownames(hyper$Sigma)
    Sigma.cnames <- colnames(hyper$Sigma)
    if(!identical(mu.names, Sigma.rnames) ||
       !identical(mu.names, Sigma.cnames)) {
      stop("names(mu), rownames(Sigma), and colnames(Sigma) are not consistent.")
    }
    if(!is.valid.vars(mu.names, c(param.names, data.names))) {
      stop("names(mu) must be a unique subset of param.names and data.names.")
    }
#    if(debug) browser()
    # indices of the variables
    var.id <- sapply(mu.names, function(x) which(x == var.names))
    # order the variables
    var.ord <- order(var.id)
    mu <- hyper$mu[var.ord]
    Sigma <- hyper$Sigma[var.ord,var.ord]
    var.id <- sort(var.id)
    # separate into theta and x components
    theta.id <- var.id[var.id <= nparams]
    x.id <- var.id[var.id > nparams] - nparams
    # format arguments
    prior.args <- list(mean = mu, cholSd = chol(Sigma),
                       thetaId = theta.id, xId = x.id)
    prior.args <- lapply(prior.args, function(x) {
      if(length(x) == 0) NULL else as.double(x)
    })
  } else {
    prior.args <- list(NULL)
  }
  prior.args
}
