#' SDE drift function.
#'
#' Computes the SDE model's drift function given data and parameter values.
#'
#' @param model An `sde.model` object.
#' @param x A vector or matrix of data with `ndims` columns.
#' @param theta A vector or matrix of parameters with `nparams` columns.
#'
#' @return A matrix with `ndims` columns containing the drift function evaluated at `x` and `theta`.  If either input contains invalid SDE data or parameters an error is thrown.
#'
#' @example examples/sde.drift.R
#' @export
sde.drift <- function(model, x, theta) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
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
#  if(debug) browser()
  # validate
  if(!all(.is.valid.params(model, theta, single.theta, nreps))) {
    stop("theta contains invalid sde parameters.")
  }
  if(!all(.is.valid.data(model, x, theta, single.x, single.theta, nreps))) {
    stop("x contains invalid sde data.")
  }
  dr <- model$cobj$Drift(xIn = as.double(x),
                         thetaIn = as.double(theta),
                         singleX = as.logical(single.x),
                         singleTheta = as.logical(single.theta),
                         nReps = as.integer(nreps))
  matrix(dr, nrow = nreps, ncol = model$ndims, byrow = TRUE)
}
