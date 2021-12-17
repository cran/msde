#' SDE prior function.
#'
#' Evaluates the SDE prior given data, parameter, and hyperparameter values.
#'
#' @param model An `sde.model` object.
#' @param theta A vector or matrix of parameters with `nparams` columns.
#' @param x A vector or matrix of data with `ndims` columns.
#' @param hyper The hyperparameters of the SDE prior.  See Details.
#'
#' @return A vector of log-prior densities evaluated at the inputs.
#'
#' @details The prior is constructed at the `C++` level by defining a function (i.e., public member)
#' ```
#' double logPrior(double *theta, double *x)
#' ```
#' within the `sdePrior` class.  At the `R` level, the `hyper.check` argument of [sde.make.model()] is a function with arguments `hyper`, `param.names`, `data.names` used to convert `hyper` into a list of `NULL` or double-vectors which get passed on to the `C++` code.  This function can also be used to throw `R`-level errors to protect the `C++` code from invalid inputs, as is done for the default prior in [mvn.hyper.check()].  For a full example see the "Custom Prior" section in `vignette("msde-quicktut")`.
#' @example examples/sde.prior.R
#' @export
sde.prior <- function(model, theta, x, hyper) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  # initialize
  param.names <- model$param.names
  data.names <- model$data.names
  x <- .format.data(x, data.names, type = "matrix")
  theta <- .format.params(theta, param.names)
  # check singles and compatible x and theta
  nreps <- c(ncol(x), ncol(theta))
  single.x <- nreps[1] == 1
  single.theta <- nreps[2] == 1
  if(!is.valid.nreps(nreps)) {
    stop("x and theta have incompatible dimensions.")
  }
  nreps <- max(nreps)
  # format hyperparameters
  phi <- model$hyper.check(hyper = hyper,
                           param.names = param.names, data.names = data.names)
  # C++ format check (is phi a list with vector-double elements)
  if(!is.valid.hyper(phi)) {
    stop("model$hyper.check must convert hyper to a list with NULL or numeric-vector elements.")
  }
  # compute
  ans <- model$cobj$Prior(thetaIn = as.double(theta),
                          xIn = as.double(x),
                          singleTheta = as.logical(single.theta),
                          singleX = as.logical(single.x),
                          nReps = as.integer(nreps),
                          phiIn = phi)
  ans
}
