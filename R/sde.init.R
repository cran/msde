#' MCMC initialization.
#'
#' Specifies the observed SDE data, interobservation times, initial parameter and missing data values to be supplied to [sde.post()].
#'
#' @param model An `sde.model` object.
#' @param x An `nobs x ndims` matrix of data.
#' @param dt A scalar or length `nobs-1` vector of interobservations times.
#' @param theta A length `nparams` vector of parameter values.
#' @param nvar.obs A scalar or length `nobs` vector of integers between 0 and `ndims` denoting the number of observed SDE variables in each row of `data`.  Defaults to `ndims`.  See Details.
#' @param m Positive integer, such that `m-1` evenly-spaced missing data time points are placed between observations.  See Details.
#'
#' @return An `sde.init` object, corresponding to a list with elements:
#' \describe{
#'   \item{`data`}{An `ncomp x ndims` matrix of complete data, where `ncomp = N_m = m * (nobs-1)+1`.}
#'   \item{`dt.m`}{The complete data interobservation time, `dt_m = dt/m`.}
#'   \item{`nvar.obs.m`}{The number of variables observed per row of `data`.  Note that `nvar.obs.m[(i-1)*m+1] == nvar.obs[ii]`, and that `nvar.obs.m[i-1] == 0` if `i` is not a multiple of `m`.}
#'   \item{`params`}{Parameter initial values.}
#' }
#' @example examples/sde.init.R
#' @export
sde.init <- function(model, x, dt, m = 1, nvar.obs, theta) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  # check inputs
  x <- .format.data(x, model$data.names, type = "matrix", strict = TRUE)
  x <- t(x)
  theta <- .format.params(theta, model$param.names)
  if(ncol(theta) > 1) {
    stop("theta must be a vector of length nparams.")
  } else {
    theta <- theta[,1]
  }
  nobs <- nrow(x)
  ndims <- model$ndims
  mm <- m
  #m <- m-1
  #if(missing(m)) m <- 2^k-1
  ncomp <- (nobs-1)*mm+1
  init.data <- matrix(NA, ncomp, ndims)
  colnames(init.data) <- model$data.names
  # interpolation to create missing data
  #if(debug) browser()
  if(length(dt) == 1) dt <- rep(dt, nobs-1)
  if(length(dt) != nobs-1) stop("x and dt have incompatible sizes.")
  if(missing(nvar.obs)) nvar.obs <- ndims
  if(length(nvar.obs) == 1) nvar.obs <- rep(nvar.obs, nobs)
  if(length(nvar.obs) != nobs) stop("nvar.obs and x have incompatible sizes.")
  if(!all(nvar.obs == round(nvar.obs)) ||
     any(nvar.obs < 0) || any(nvar.obs > ndims)) {
    stop("nvar.obs must be an integer scalar or vector with values between 0 and ndims.")
  }
  dtnew <- rep(dt/mm, each = mm)
  told <- cumsum(c(0, dt))
  tnew <- cumsum(c(0, dtnew))
  for(ii in 1:ndims) {
    init.data[,ii] <- approx(x = told, y = x[,ii],
                             xout = tnew)$y
  }
  # validate
  if(!all(.is.valid.params(model, t(theta), TRUE, 1))) {
    stop("Invalid sde parameters.")
  }
  if(!all(.is.valid.data(model, t(init.data), t(theta), FALSE, TRUE, ncomp))) {
    stop("Invalid sde complete data.")
  }
  #if(dt1) dtnew <- dtnew[1]
  #par.index <- nvar.obs
  #if(missing(par.index)) par.index <- ndims
  par.ind <- rep(0, ncomp)
  par.ind[seq(1, ncomp, len = nobs)] <- nvar.obs
  ans <- list(data = init.data, dt.m = dtnew,
              nvar.obs.m = par.ind, params = theta)
  #if(!missing(params)) ans$params <- params
  class(ans) <- "sde.init"
  ans
}
