#' SDE loglikelihood function.
#'
#' Evaluates the loglikelihood function given SDE data and parameter values.
#' @param model An \code{sde.model} object.
#' @param x A matrix or 3-d array of data with \code{dim(x)[1]} observations and \code{dim(x)[2] == ndims}.
#' @param dt A scalar or vector of length \code{dim(x)[1]-1} of time intervals between observations.
#' @param theta A vector or matrix of parameters with \code{nparams} columns.
#' @param ncores If \code{model} is compiled with \code{OpenMP}, the number of cores to use for parallel processing.  Otherwise, uses \code{ncores = 1} and gives a warning.
#' @return A vector of loglikelihood evaluations, of the same length as the  third dimension of \code{x} and/or first dimension of \code{theta}.  If input contains invalid data or parameters an error is thrown.
#' @examples
#' \donttest{
#' # Compile Heston's model
#' hex <- example.models("hest")
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#'
#' # Simulate data
#' nreps <- 10
#' nobs <- 100
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' Theta <- apply(t(replicate(nreps, theta)), 2, jitter)
#' x0 <- c(X = log(1000), Z = 0.1)
#' X0 <- apply(t(replicate(nreps,x0)), 2, jitter)
#' dT <- 1/252
#' hsim <- sde.sim(model = hmod, x0 = X0, theta = Theta,
#'                 dt = dT, dt.sim = dT/10, nobs = nobs, nreps = nreps)
#'
#' # single parameter, single data
#' sde.loglik(model = hmod, x = hsim$data[,,1], dt = dT, theta = theta)
#' # multiple parameters, single data
#' sde.loglik(model = hmod, x = hsim$data[,,1], dt = dT, theta = Theta)
#' # multiple parameters, multiple data
#' sde.loglik(model = hmod, x = hsim$data, dt = dT, theta = Theta)
#' }
#' @export
sde.loglik <- function(model, x, dt, theta, ncores = 1) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  x <- .format.data(x, model$data.names, type = "array")
  theta <- .format.params(theta, model$param.names)
  # problem dimensions
  ncomp <- dim(x)[2]
  if(ncomp <= 2) {
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
  model$loglik(xIn = as.double(x), dTIn = as.double(dt),
               thetaIn = as.double(theta),
               nComp = as.integer(ncomp),
               nReps = as.integer(nreps),
               singleX = as.logical(single.x),
               singleTheta = as.logical(single.theta),
               nCores = as.integer(ncores))
}
