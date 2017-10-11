#' SDE data and parameter validators.
#'
#' Checks whether input SDE data and parameters are valid.
#' @param model An \code{sde.model} object.
#' @param x A length-\code{ndims} vector or \code{ndims}-column matrix of SDE data.
#' @param theta A length-\code{nparams} vector or \code{nparams}-column of SDE parameter values.
#' @return A logical scalar or vector indicating whether the given data/parameter pair is valid.
#' @name sde.valid
#' @examples
#' # Heston's model
#' # valid data is: Z > 0
#' # valid parameters are: gamma, sigma > 0, |rho| < 1, beta > .5 * sigma^2
#' hmod <- sde.examples("hest") # load model
#'
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#'
#' # valid data
#' x0 <- c(X = log(1000), Z = 0.1)
#' sde.valid.data(model = hmod, x = x0, theta = theta)
#'
#' # invalid data
#' x0 <- c(X = log(1000), Z = -0.1)
#' sde.valid.data(model = hmod, x = x0, theta = theta)
#'
#' # valid parameters
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' sde.valid.params(model = hmod, theta = theta)
#'
#' # invalid parameters
#' theta <- c(alpha = 0.1, gamma = -4, beta = 0.8, sigma = 0.6, rho = -0.8)
#' sde.valid.params(model = hmod, theta = theta)
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
  .sde_isData(sdeptr = model$ptr, xIn = as.double(x),
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
  .sde_isParams(sdeptr = model$ptr, thetaIn = as.double(theta),
                nReps = as.integer(nreps))
}

#--- internal versions: no argument checking/formatting -------------------

.is.valid.data <- function(model, x, theta, single.x, single.theta,
                           nreps) {
  .sde_isData(sdeptr = model$ptr, xIn = as.double(x),
              thetaIn = as.double(theta),
              singleX = as.logical(single.x),
              singleTheta = as.logical(single.theta),
              nReps = as.integer(nreps))
}

.is.valid.params <- function(model, theta, single.theta, nreps) {
  rep(.sde_isParams(sdeptr = model$ptr, thetaIn = as.double(theta),
                    nReps = as.integer(ifelse(single.theta, 1, nreps))),
      times = ifelse(single.theta, nreps, 1))
}
