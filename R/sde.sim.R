#' Simulation of multivariate SDE trajectories.
#'
#' Simulates a discretized Euler-Maruyama approximation to the true SDE trajectory.
#' @param model An `sde.model` object.
#' @param x0 A vector or a matrix of size `nreps x ndims` of the SDE values at time 0.
#' @param theta A vector or matrix of size `nreps x nparams` of SDE parameters.
#' @param dt Scalar interobservation time.
#' @param dt.sim Scalar interobservation time for simulation.  That is, interally the interobservation time is `dt.sim` but only one out of every `dt/dt.sim` simulation steps is kept in the output.
#' @param nobs The number of SDE observations per trajectory to generate.
#' @param burn Scalar burn-in value.  Either an integer giving the number of burn-in steps, or a value between 0 and 1 giving the fraction of burn-in relative to `nobs`.
#' @param nreps The number of SDE trajectories to generate.
#' @param max.bad.draws The maximum number of times that invalid forward steps are proposed.  See Details.
#' @param verbose Whether or not to display information on the simulation.
#' @details The simulation algorithm is a Markov process with \eqn{Y_0 = x_0} and
#' \deqn{
#' Y_{t+1} \sim \mathcal{N}(Y_t + \mathrm{dr}(Y_t, \theta) dt_{\mathrm{sim}}, \mathrm{df}(Y_t, \theta) dt_{\mathrm{sim}}),
#' }{
#' Y_(t+1) ~ N(Y_t + dr(Y_t, \theta) dt_(sim), df(Y_t, \theta) dt_(sim)),
#'}
#' where \eqn{\mathrm{dr}(y, \theta)}{dr(y, \theta)} is the SDE drift function and \eqn{\mathrm{df}(y, \theta)}{df(y, \theta)} is the diffusion function on the **variance** scale.  At each step, a while-loop is used until a valid SDE draw is produced.  The simulation algorithm terminates after `nreps` trajectories are drawn or once a total of `max.bad.draws` are reached.
#' @return A list with elements:
#' \describe{
#'   \item{`data`}{An array of size `nobs x ndims x nreps` containing the simulated SDE trajectories.}
#'   \item{`params`}{The vector or matrix of parameter values used to generate the data.}
#'   \item{`dt, dt.sim`}{The actual and internal interobservation times.}
#'   \item{`nbad`}{The total number of bad draws.}
#' }
#' @example examples/sde.sim.R
#' @export
sde.sim <- function(model, x0, theta, dt, dt.sim,
                    nobs, burn = 0, nreps = 1,
                    max.bad.draws = 5e3, verbose = TRUE) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  # format data and parameters
  ndims <- model$ndims
  x0 <- .format.data(x0, model$data.names, type = "matrix")
  theta <- .format.params(theta, model$param.names)
  # check validity
  N <- c(ncol(x0), ncol(theta))
  single.x <- N[1] == 1
  single.theta <- N[2] == 1
  if(!is.valid.nreps(c(N, nreps))) {
    stop("Incompatible dimensions of x0, theta, and nreps.")
  }
  if(!all(.is.valid.params(model, theta, single.theta, nreps))) {
    stop("theta contains invalid sde parameters.")
  }
  if(!all(.is.valid.data(model, x0, theta, single.x, single.theta, nreps))) {
    stop("x0 contains invalid sde data.")
  }
  # time
  if(dt.sim <= dt) {
    rr <- ceiling(dt/dt.sim)
    dT <- dt/rr
  } else {
    rr <- 1
    dT <- dt
  }
  if(burn < 1) burn <- nobs*burn
  burn <- floor(burn)
  if(verbose) {
    message("Number of normal draws required: ",
            round((nobs+burn)*rr*nreps, 2))
    message("Running simulation...")
  }
#  if(debug) browser()
  tm <- chrono()
  ans <- model$cobj$Sim(nDataOut = as.integer(nobs*ndims*nreps),
                        N = as.integer(nobs),
                        burn = as.integer(burn),
                        reps = as.integer(nreps),
                        r = as.integer(rr),
                        dT = as.double(dT),
                        MAXBAD = as.integer(max.bad.draws),
                        initData = as.double(x0),
                        params = as.double(theta),
                        singleX = as.logical(single.x),
                        singleTheta = as.logical(single.theta))
  tm <- chrono(tm, display = verbose)
  names(ans) <- c("dataOut", "nBadDraws")
  if(verbose) message("Bad Draws = ", ans$nBadDraws)
  data <- aperm(array(ans$dataOut, dim = c(ndims, nobs, nreps)),
                perm = c(2,1,3))
  dimnames(data) <- list(NULL, model$data.names, NULL)
  out <- list(data = drop(data), params = drop(t(theta)),
              dt = dt, dt.sim = dT,
              nbad = ans$nBadDraws)
  out
}
