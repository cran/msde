#' MCMC initialization.
#'
#' Specifies the observed SDE data, interobservation times, initial parameter and missing data values to be supplied to \code{\link{sde.post}}.
#' @param model An \code{sde.model} object.
#' @param x An \code{nobs x ndims} matrix of data.
#' @param dt A scalar or length \code{nobs-1} vector of interobservations times.
#' @param theta A length \code{nparams} vector of parameter values.
#' @param nvar.obs A scalar or length \code{nobs} vector of integers between 0 and \code{ndims} denoting the number of observed SDE variables in each row of \code{data}.  Defaults to \code{ndims}.  See Details.
#' @param m Positive integer, such that \code{m-1} evenly-spaced missing data time points are placed between observations.  See Details.
#' @return An \code{sde.init} object, corresponding to a list with elements:
#' \describe{
#'   \item{\code{data}}{An \code{ncomp x ndims} matrix of complete data, where \code{ncomp = N_m = m * (nobs-1)+1}.}
#'   \item{\code{dt.m}}{The complete data interobservation time, \code{dt_m = dt/m}.}
#'   \item{\code{nvar.obs.m}}{The number of variables observed per row of \code{data}.  Note that \code{nvar.obs.m[(i-1)*m+1] == nvar.obs[ii]}, and that \code{nvar.obs.m[i-1] == 0} if \code{i} is not a multiple of \code{m}.}
#'   \item{\code{params}}{Parameter initial values.}
#' }
#' @importFrom stats approx
#' @examples
#' # load Heston's model
#' hmod <- sde.examples("hest")
#'
#' # generate some observed data
#' nObs <- 5
#' x0 <- c(X = log(1000), Z = 0.1)
#' X0 <- apply(t(replicate(nObs, x0)), 2, jitter)
#' dT <- .6
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#'
#' # no missing data
#' sde.init(model = hmod, x = X0, dt = dT, theta = theta)
#'
#' # all but endpoint volatilities are missing
#' sde.init(model = hmod, x = X0, dt = dT, m = 1,
#'          nvar.obs = c(2, rep(1, nObs-2), 2), theta = theta)
#'
#' # all volatilities missing,
#' # two completely missing SDE timepoints between observations
#' m <- 3 # divide each observation interval into m equally spaced timepoints
#' sde.init(model = hmod, x = X0, dt = dT,
#'          m = m, nvar.obs = 1, theta = theta)
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
