#' SDE loglikelihood function.
#'
#' Evaluates the loglikelihood function given SDE data and parameter values.
#'
#' @param model An `sde.model` object.
#' @param x A matrix or 3-d array of data with `dim(x)[1]` observations and `dim(x)[2] == ndims`.
#' @param dt A scalar or vector of length `dim(x)[1]-1` of time intervals between observations.
#' @param theta A vector or matrix of parameters with `nparams` columns.
#' @param ncores If `model` is compiled with `OpenMP`, the number of cores to use for parallel processing.  Otherwise, uses `ncores = 1` and gives a warning.
#'
#' @return A vector of loglikelihood evaluations, of the same length as the  third dimension of `x` and/or first dimension of `theta`.  If input contains invalid data or parameters an error is thrown.
#'
#' @example examples/sde.loglik.R
#' @export
sde.loglik <- function(model, x, dt, theta, ncores = 1) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  x <- .format.data(x, model$data.names, type = "array")
  theta <- .format.params(theta, model$param.names)
  # problem dimensions
  ncomp <- dim(x)[2]
  if(ncomp < 2) {
    stop("likelihood calculation requires at least two observations.")
  }
  if(length(dt) == 1) dt <- rep(dt, ncomp-1)
  if(length(dt) != ncomp-1) {
    stop("x and dt have incompatible dimensions.")
  }
  nreps <- c(dim(x)[3], ncol(theta))
  single.x <- nreps[1] == 1
  single.theta <- nreps[2] == 1
  if(!is.valid.nreps(nreps)) {
    stop("x and theta have incompatible dimensions.")
  }
  nreps <- max(nreps)
  # multicore functionality
  if(ncores < 1) stop("ncores must be a positive integer.")
  if(!model$omp && ncores > 1) {
    warning("model not compiled with openMP: ncores set to 1.")
    ncores <- 1
  }
  # validate
  if(!all(.is.valid.params(model, theta, single.theta, nreps))) {
    stop("theta contains invalid sde parameters.")
  }
  if(!single.theta) {
    theta.ind <- rep(1:nreps, each = ncomp)
  } else {
    theta.ind <- 1
  }
  if(!all(.is.valid.data(model, x, theta[,theta.ind],
                         single.x, single.theta, nreps*ncomp))) {
    stop("x contains invalid sde data.")
  }
  # compute
  model$cobj$Loglik(xIn = as.double(x),
                    dTIn = as.double(dt),
                    thetaIn = as.double(theta),
                    nComp = as.integer(ncomp),
                    nReps = as.integer(nreps),
                    singleX = as.logical(single.x),
                    singleTheta = as.logical(single.theta),
                    nCores = as.integer(ncores))
}
