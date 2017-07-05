#' MCMC sampler for the SDE posterior.
#'
#' A Metropolis-within-Gibbs sampler for the Euler-Maruyama approximation to the true posterior density.
#' @param model An \code{sde.model} object constructed with \code{\link{sde.make.model}}.
#' @param init An \code{sde.init} object constructed with \code{\link{sde.init}}.
#' @param hyper The hyperparameters of the SDE prior.  See \code{\link{sde.prior}}.
#' @param nsamples Number of MCMC iterations.
#' @param burn Integer number of burn-in samples, or fraction of \code{nsamples} to prepend as burn-in.
#' @param mwg.sd Standard deviation jump size for Metropolis-within-Gibbs on parameters and missing components of first SDE observation (see Details).
#' @param adapt Logical or list to specify adaptive Metropolis-within-Gibbs sampling (see Details).
#' @param loglik.out Logical, whether to return the loglikelihood at each step.
#' @param last.miss.out Logical, whether to return the missing sde components of the last observation.
#' @param update.data Logical, whether to update the missing data.
#' @param data.out A scalar, integer vector, or list of three integer vectors determining the subset of data to be returned (see Details).
#' @param update.params Logical, whether to update the model parameters.
#' @param fixed.params Logical vector of length \code{nparams} indicating which parameters are to be held fixed in the MCMC sampler.
#' @param ncores If \code{model} is compiled with \code{OpenMP}, the number of cores to use for parallel processing.  Otherwise, uses \code{ncores = 1} and gives a warning.
#' @param verbose Logical, whether to periodically output MCMC status.
#' @details The Metropolis-within-Gibbs (MWG) jump sizes can be specified as a scalar, a vector or length \code{nparams + ndims}, or a named vector containing the elements defined by \code{sde.init$nvar.obs.m[1]} (the missing variables in the first SDE observation) and \code{fixed.params} (the SDE parameters which are not held fixed).  The default jump sizes for each MWG random variable are \code{.25 * |initial_value|}.
#'
#' \code{adapt == TRUE} implements an adaptive MCMC proposal by Rosenthal and Roberts (2005).  At step \eqn{n} of the MCMC, the jump size of each MWG random variable is increased or decreased by \eqn{\delta(n)}, depending on whether the cumulative acceptance rate is above or below the optimal value of 0.44.  If \eqn{\sigma_n} is the size of the jump at step \eqn{n}, then the next jump size is determined by
#' \deqn{
#' \log(\sigma_{n+1}) = \log(\sigma_n) \pm \delta(n), \qquad \delta(n) = \min(.01, 1/n^{1/2}).
#' }{
#' log(sigma_(n+1)) = log(sigma_n) \pm delta(n), where delta(n) = min(.01, 1/n^(1/2)).
#' }
#' When \code{adapt} is not logical, it is a list with elements \code{max} and \code{rate}, such that \code{delta(n) = min(max, 1/n^rate)}.  These elements can be scalars or vectors in the same manner as \code{mwg.sd}.
#'
#' For SDE models with thousands of latent variables, \code{data.out} can be used to thin the MCMC missing data output.  An integer vector or scalar returns specific or evenly-spaced posterior samples from the \code{ncomp x ndims} complete data matrix.  A list with elements \code{isamples}, \code{icomp}, and \code{idims} determines which samples, time points, and SDE variables to return.  The first of these can be a scalar or vector with the same meaning as before.
#' @return A list with elements:
#' \describe{
#'   \item{\code{params}}{An \code{nsamples x nparams} matrix of posterior parameter draws.}
#'   \item{\code{data}}{A 3-d array of posterior missing data draws, for which the output dimensions are specified by \code{data.out}.}
#'   \item{\code{data.out}}{A list of three integer vectors specifying which timepoints, variables, and MCMC iterations correspond to the values in the \code{data} output.}
#'   \item{\code{init}}{The \code{sde.init} object which initialized the sampler.}
#'   \item{\code{loglik}}{If \code{loglik.out == TRUE}, the vector of \code{nsamples} complete data loglikelihoods calculated at each posterior sample.}
#'   \item{\code{last.iter}}{A list with elements \code{data} and \code{params} giving the last MCMC sample.  Useful for resuming the MCMC from that point.}
#'   \item{\code{last.miss}}{If \code{last.miss.out == TRUE}, an \code{nsamples x nmissN} matrix of all posterior draws for the missing data in the final observation.  Useful for SDE forecasting at future timepoints.}
#'   \item{\code{accept}}{A named list of acceptance rates for the various components of the MCMC sampler.}
#' }
#' @examples
#' \donttest{
#' # Posterior inference for Heston's model
#' hex <- example.models("hest")
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#'
#' # Simulate data
#' X0 <- c(X = log(1000), Z = 0.1)
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' dT <- 1/252
#' nobs <- 1000
#' hest.sim <- sde.sim(model = hmod, x0 = X0, theta = theta,
#'                     dt = dT, dt.sim = dT/10, nobs = nobs)
#'
#' # initialize MCMC sampler
#' # both components observed, no missing data between observations
#' init <- sde.init(model = hmod, x = hest.sim$data,
#'                  dt = hest.sim$dt, theta = theta)
#'
#' # Initialize posterior sampling argument
#' nsamples <- 1e4
#' burn <- 1e3
#' hyper <- NULL # flat prior
#' hest.post <- sde.post(model = hmod, init = init, hyper = hyper,
#'                       nsamples = nsamples, burn = burn)
#'
#' # plot the histogram for the sampled parameters
#' par(mfrow = c(2,3))
#' for(ii in 1:length(hmod$param.names)) {
#'   hist(hest.post$params[,ii],breaks=100, freq = FALSE,
#'        main = parse(text = hmod$param.names[ii]), xlab = "")
#' }
#' }
#' @export
sde.post <- function(model, init, hyper,
                     nsamples, burn, mwg.sd = NULL, adapt = TRUE,
                     loglik.out = FALSE, last.miss.out = FALSE,
                     update.data = TRUE, data.out,
                     update.params = TRUE, fixed.params,
                     ncores = 1, verbose = TRUE) {
  # model constants
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initial values
  if(class(init) != "sde.init") {
    # all argument checking done at construction time
    stop("init must be an sde.init object.")
  }
  dt <- init$dt.m
  par.index <- init$nvar.obs.m
  init.data <- init$data
  init.params <- init$params
  if(missing(fixed.params)) fixed.params <- rep(FALSE, nparams)
  # parse inputs
  ncomp <- nrow(init.data)
  nmiss0 <- ndims - par.index[1]
  nparams2 <- nparams+nmiss0
  ## if(length(dt) == 1) dt <- rep(dt, ncomp-1) # time
  ## .check.init(init.data, dt, init.params, param.names, data.names)
  if(missing(burn)) burn <- max(.1, 1e3)
  if(burn < 1) burn <- nsamples*burn
  burn <- floor(burn)
  # output sizes
  if(all(par.index == ndims)) update.data <- FALSE
  data.out <- .set.data.out(data.out, nsamples, ncomp, ndims, update.data)
  nsamples.out <- length(data.out$isamples)
  ncomp.out <- length(data.out$icomp)
  ndims.out <- length(data.out$idims)
  ndata.out <- ifelse(update.data, nsamples.out*ndims.out*ncomp.out, 1)
  nparams.out <- ifelse(update.params, nsamples*nparams, 1)
  nmissN <- ndims - par.index[ncomp] # last missing output
  if(nmissN == 0) last.miss.out <- FALSE
  nlast.miss.out <- ifelse(last.miss.out, nsamples*nmissN, 1)
  # prior specification
  prior <- hyper
  # format hyperparameters
  prior <- model$hyper.check(prior, param.names, data.names)
  # C++ format check (is phi a list with vector-double elements)
  if(!is.valid.hyper(prior)) {
    stop("model$hyper.check must convert hyper to a list with NULL or vector-numeric elements.")
  }
  # random walk jump size
  if(is.null(mwg.sd)) mwg.sd <- .25 * abs(c(init.params, init.data[1,]))
  tune.par <- .set.jump(mwg.sd, adapt, param.names, data.names)
  # multicore functionality
  if(ncores < 1) stop("ncores must be a positive integer.")
  if(!model$omp && ncores > 1) {
    warning("model not compiled with openMP: ncores set to 1.")
    ncores <- 1
  }
  if(verbose) {
    message("Output size:")
    if(update.params) message("params = ", round(nparams.out, 2))
    if(update.data) message("data = ", round(ndata.out, 2))
    message("Running posterior sampler...")
  }
#  if(debug) browser()
  tm <- chrono()
  # compute
  ans <- model$post(initParams = as.double(init.params),
                    initData = as.double(t(init.data)),
                    dT = as.double(dt),
                    nDimsPerObs = as.integer(par.index),
                    fixedParams = as.logical(fixed.params),
                    nSamples = as.integer(nsamples),
                    burn = as.integer(burn),
                    nParamsOut = as.integer(nparams.out),
                    nDataOut = as.integer(ndata.out),
                    dataOutSmp = as.integer(data.out$isamples-1),
                    dataOutComp = as.integer(data.out$icomp-1),
                    dataOutDims = as.integer(data.out$idims-1),
                    updateParams = as.double(update.params),
                    updateData = as.double(update.data),
                    priorArgs = prior,
                    tunePar = tune.par,
                    #rwJumpSd = as.double(tune.par$sd),
                    updateLogLik = as.integer(loglik.out),
                    nLogLikOut = as.integer(ifelse(loglik.out, nsamples, 1)),
                    updateLastMiss = as.integer(last.miss.out),
                    nLastMissOut = as.integer(nlast.miss.out),
                    nCores = as.integer(ncores))
  tm <- chrono(tm, display = verbose)
  names(ans) <- c("paramsOut", "dataOut", "paramAccept", "gibbsAccept",
                  "logLikOut", "lastMissOut", "lastIter", "mwgSd")
  # acceptance rates
#  if(debug) browser()
  accept <- .set.accept(ans$gibbsAccept, ans$paramAccept,
                        nsamples+burn, par.index, fixed.params,
                        param.names, data.names,
                        update.data, update.params, verbose)
  out <- list()
  if(update.params) {
    theta.out <- matrix(ans$paramsOut,
                        ncol = nparams, byrow = TRUE,
                        dimnames = list(NULL, param.names))
    out <- c(out, list(params = theta.out))
  } else out <- c(out, list(params = init.params))
  if(update.data) {
    x.out <- array(ans$dataOut,
                   dim = c(ndims.out, ncomp.out, nsamples.out),
                   dimnames = list(data.names[data.out$idims], NULL, NULL))
    x.out <- aperm(x.out, perm = c(2,1,3))
    out <- c(out, list(data = x.out))
  } else out <- c(out, list(data = init.data))
  if(loglik.out) out <- c(out, list(loglik = ans$logLikOut))
  out <- c(out, list(init = init, data.out = data.out,
                     mwg.sd = ans$mwgSd, hyper = hyper))
  last.iter <- list(params = ans$lastIter[1:nparams],
                    data = matrix(ans$lastIter[nparams + 1:(ncomp*ndims)],
                                  ncomp, ndims, byrow = TRUE))
  names(last.iter$params) <- param.names
  colnames(last.iter$data) <- data.names
  out <- c(out, list(last.iter = last.iter))
  if(last.miss.out) {
    last.miss <- matrix(ans$lastMissOut, nsamples, nmissN, byrow = TRUE)
    colnames(last.miss) <- data.names[par.index[ncomp]+(1:nmissN)]
    out <- c(out, list(last.miss = last.miss))
  }
  out <- c(out, list(accept = accept))
  out
}

#--- helper functions ----------------------------------------------------------

# which MCMC iterations and which time points are returned
# input is a scalar, integer vector, or list of integer vectors named
# with names(data.out) = c("icomp", "idims", "isamples") (in any order)
# output is a list of integer vectors with elements named as above
.set.data.out <- function(data.out, nsamples, ncomp, ndims, update.data) {
  if(missing(data.out)) data.out <- 2e3
  if(!is.list(data.out)) {
    # keep complete data matrix at isamples == data.out
    isamples <- data.out
    icomp <- 1:ncomp
    idims <- 1:ndims
  } else {
    if(!identical(sort(names(data.out)),
                  sort(c("icomp", "idims", "isamples")))) {
      stop("data.out must be scalar, vector, or list with elements icomp, idims, and isamples.")
    }
    isamples <- data.out$isamples
    icomp <- data.out$icomp
    idims <- data.out$idims
  }
  if(length(isamples) == 1) {
    # evenly space returned samples
   isamples <- unique(floor(seq(1, nsamples, len = isamples)))
  }
  if(!update.data) {
    # return the data once
    isamples <- 1:nsamples
    icomp <- 1:ncomp
    idims <- 1:ndims
  }
  # check inputs
  if(anyDuplicated(icomp) || !all(icomp %in% 1:ncomp)) {
    stop("data.out$icomp must be a vector of integers between 1 and ncomp.")
  }
  if(anyDuplicated(idims) || !all(idims %in% 1:ndims)) {
    stop("data.out$idims must be a vector of integers between 1 and ndims.")
  }
  if(anyDuplicated(isamples) || !all(isamples %in% 1:nsamples)) {
    stop("data.out$isamples must be a vector of integers between 1 and nsamples.")
  }
  list(icomp = sort(icomp), idims = sort(idims), isamples = sort(isamples))
}

# jump sizes
.set.jump <- function(mwg.sd, adapt, param.names, data.names) {
  nparams <- length(param.names)
  ndims <- length(data.names)
  .format.arg <- function(x) {
    if(length(x) == 1) {
      x <- rep(x, nparams+ndims)
      names(x) <- c(param.names, data.names)
    } else if(length(x) == nparams+ndims && is.null(names(x))) {
      names(x) <- c(param.names, data.names)
    }
    x
  }
  .set.arg <- function(x, id) {
    y <- rep(0, nparams+ndims)
    names(y) <- c(param.names, data.names)
    y[names(x)] <- x
    y
  }
  mwg.sd <- .format.arg(mwg.sd)
  if(!is.valid.vars(names(mwg.sd), c(param.names, data.names))) {
    stop("names(mwg.sd) must be a unique subset of param.names and data.names.")
  }
  mwg.sd <- .set.arg(mwg.sd)
  # adaptive MCMC
  if(is.logical(adapt)) {
    if(adapt) {
      amax <- .01
      arate <- .5
    } else {
      amax <- 0
      arate <- 0
    }
  } else {
    adapt <- TRUE
    amax <- adapt$max
    arate <- adapt$rate
    if(is.null(amax) || is.null(arate)) {
      stop("adapt must be logical or a list with elements max and rate.")
    }
  }
  # check args
  amax <- .format.arg(amax)
  if(!is.valid.vars(names(amax), c(param.names, data.names))) {
    stop("names(adapt$max) must be a unique subset of param.names and data.names.")
  }
  amax <- .set.arg(amax)
  arate <- .format.arg(arate)
  if(!is.valid.vars(names(arate), c(param.names, data.names))) {
    stop("names(adapt$rate) must be a unique subset of param.names and data.names.")
  }
  arate <- .set.arg(arate)
  if(any(arate < 0)) stop("adapt$rate must be non-negative.")
  list(sd = as.double(mwg.sd), max = as.double(amax),
       rate = as.double(arate), adapt = as.logical(amax > 0))
}

# parse acceptance rates
.set.accept <- function(bb.accept, vnl.accept, nsamples,
                        par.index, fixed.params,
                        param.names, data.names,
                        update.data, update.params, verbose) {
  ndims <- length(data.names)
  nparams <- length(param.names)
  accept <- NULL
  if(update.data) {
    accept <- c(accept, list(data = bb.accept[-1]/nsamples))
    bb.acc <- accept$data[par.index[-1] < ndims]*100
    bb.acc <- signif(c(min(bb.acc), mean(bb.acc)), 3)
    if(verbose) {
      message("Bridge accept: min = ", bb.acc[1],
              "%, avg = ", bb.acc[2], "%")
    }
  }
  if(update.params) {
    accept <- c(accept, list(params = vnl.accept[1:nparams]/nsamples))
    if(verbose) {
      for(ii in which(!fixed.params)) {
        message(param.names[ii], " accept: ",
                signif(accept$params[ii]*100,3), "%")
      }
    }
  }
  nmiss0 <- ndims-par.index[1]
  if(update.data && (nmiss0 > 0)) {
    accept <- c(accept,
                list(miss0 = vnl.accept[nparams+(1:nmiss0)]/nsamples))
    for(ii in 1:nmiss0) {
      message(data.names[par.index[1] + ii], "0 accept: ",
              signif(accept$miss0[ii]*100,3), "%")
    }
  }
  accept
}
