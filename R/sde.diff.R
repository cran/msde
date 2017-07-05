#' SDE diffusion function.
#'
#' Computes the SDE model's diffusion function given data and parameter values.
#' @param model An \code{sde.model} object.
#' @param x A vector or matrix of data with \code{ndims} columns.
#' @param theta A vector or matrix of parameters with \code{nparams} columns.
#' @return A matrix with \code{ndims^2} columns containing the diffusion function evaluated at \code{x} and \code{theta}. Each row corresponds to the upper triangular Cholesky factor of the diffusion matrix.  If either input contains invalid SDE data or parameters an error is thrown.
#' @examples
#' \donttest{
#' # compile model
#' hex <- example.models("hest")
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#'
#' # single input
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' x0 <- c(X = log(1000), Z = 0.1)
#' sde.diff(model = hmod, x = x0, theta = theta)
#'
#' # multiple inputs
#' nreps <- 10
#' Theta <- apply(t(replicate(nreps, theta)), 2, jitter)
#' X0 <- apply(t(replicate(nreps, x0)), 2, jitter)
#' sde.diff(model = hmod, x = X0, theta = Theta)
#'
#' # mixed inputs
#' sde.diff(model = hmod, x = x0, theta = Theta)
#' }
#' @export
sde.diff <- function(model, x, theta) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  # initialize
  ndims <- model$ndims
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
  # validate
  if(!all(.is.valid.params(model, theta, single.theta, nreps))) {
    stop("theta contains invalid sde parameters.")
  }
  if(!all(.is.valid.data(model, x, theta, single.x, single.theta, nreps))) {
    stop("x contains invalid sde data.")
  }
  df <- model$diff(xIn = as.double(x),
                   thetaIn = as.double(theta),
                   singleX = as.logical(single.x),
                   singleTheta = as.logical(single.theta),
                   nReps = as.integer(nreps))
  df <- matrix(df, nrow = nreps, ncol = ndims^2, byrow = TRUE)
  # put zeros into the off-triangular elements
  df[,lower.tri(diag(ndims))] <- 0
  df
}