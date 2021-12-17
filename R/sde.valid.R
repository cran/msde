#' SDE data and parameter validators.
#'
#' Checks whether input SDE data and parameters are valid.
#' @param model An `sde.model` object.
#' @param x A length-`ndims` vector or `ndims`-column matrix of SDE data.
#' @param theta A length-`nparams` vector or `nparams`-column of SDE parameter values.
#' @return A logical scalar or vector indicating whether the given data/parameter pair is valid.
#' @name sde.valid
#' @aliases sde.valid.data sde.valid.params
#' @example examples/sde.valid.R
#' @export
sde.valid.data <- function(model, x, theta) {
  if(class(model) != "sde.model")
    stop("model must be an sde.model object.")
  # initialize
  x <- .format.data(x, model$data.names, type = "matrix")
  theta <- .format.params(theta, model$param.names)
  # check singles and compatible x and theta
  nreps <- c(ncol(x), ncol(theta))
  single.x <- nreps[1] == 1
  single.theta <- nreps[2] == 1
  if(!is.valid.nreps(nreps)) {
    stop("x and theta have incompatible dimensions.")
  }
  nreps <- max(nreps)
  model$cobj$isData(xIn = as.double(x),
                    thetaIn = as.double(theta),
                    singleX = as.logical(single.x),
                    singleTheta = as.logical(single.theta),
                    nReps = as.integer(nreps))
}

#' @rdname sde.valid
#' @export
sde.valid.params <- function(model, theta) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # initialize
  theta <- .format.params(theta, model$param.names)
  # check singles and compatible x and theta
  nreps <- ncol(theta)
  model$cobj$isParams(thetaIn = as.double(theta),
                      nReps = as.integer(nreps))
}

#--- internal versions: no argument checking/formatting -------------------

.is.valid.data <- function(model, x, theta, single.x, single.theta,
                           nreps) {
  model$cobj$isData(xIn = as.double(x),
                    thetaIn = as.double(theta),
                    singleX = as.logical(single.x),
                    singleTheta = as.logical(single.theta),
                    nReps = as.integer(nreps))
}

.is.valid.params <- function(model, theta, single.theta, nreps) {
  rep(model$cobj$isParams(thetaIn = as.double(theta),
                          nReps = as.integer(ifelse(single.theta, 1, nreps))),
      times = ifelse(single.theta, nreps, 1))
}
