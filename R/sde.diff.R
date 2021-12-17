#' SDE diffusion function.
#'
#' Computes the SDE model's diffusion function given data and parameter values.
#'
#' @param model An `sde.model` object.
#' @param x A vector or matrix of data with `ndims` columns.
#' @param theta A vector or matrix of parameters with `nparams` columns.
#'
#' @return A matrix with `ndims^2` columns containing the diffusion function evaluated at `x` and `theta`. Each row corresponds to the upper triangular Cholesky factor of the diffusion matrix.  If either input contains invalid SDE data or parameters an error is thrown.
#'
#' @example examples/sde.diff.R
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
  df <- model$cobj$Diff(xIn = as.double(x),
                        thetaIn = as.double(theta),
                        singleX = as.logical(single.x),
                        singleTheta = as.logical(single.theta),
                        nReps = as.integer(nreps))
  df <- matrix(df, nrow = nreps, ncol = ndims^2, byrow = TRUE)
  # put zeros into the off-triangular elements
  df[,lower.tri(diag(ndims))] <- 0
  df
}
