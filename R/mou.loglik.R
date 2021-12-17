#' Loglikelihood for multivariate Ornstein-Uhlenbeck process.
#'
#' Computes the exact Euler loglikelihood for any amount of missing data using a Kalman filter.
#'
#' @param X An `nobs x ndims` matrix of complete data.
#' @param dt A scalar or length `nobs-1` vector of interobservations times.
#' @param nvar.obs A scalar or length `nobs` vector of integers between 0 and `ndims` denoting the number of observed SDE variables in each row of `data`.  Defaults to `ndims`.  See [sde.init()] for details.
#' @param Gamma A `ndims x ndims` of linear-drift parameters.  See Details.
#' @param Lambda A length-`ndims` vector of constant-drift parameters.  See Details.
#' @param Phi A `ndims x ndims` positive definite variance matrix.  See Details.
#' @param mu0,Sigma0 Mean and variance of marginal multivariate normal distribution of `X[1,]`.  Defaults to iid standard normals for each component.
#'
#' @return Scalar value of the loglikelihood.  See Details.
#'
#' @details The \eqn{p}-dimensional multivariate Ornstein-Uhlenbeck (mOU) process \eqn{Y_t = (Y_{1t}, \ldots, Y_{dt})}{Y_t = (Y_1t, \ldots, Y_dt)} satisfies the SDE
#' \deqn{
#' dY_t = (\Gamma Y_t + \Lambda)dt + \Phi^{1/2} dB_t,
#' }{
#' dY_t = (\Gamma Y_t + \Lambda)dt + \Phi^(1/2) dB_t,
#' }
#' where \eqn{B_t = (B_{1t}, \ldots, B_{pt})}{B_t = (B_1t, \ldots, B_pt)} is \eqn{p}-dimensional Brownian motion.  Its Euler discretization is of the form
#' \deqn{
#' Y_{n+1} = Y_n + (\Gamma Y_n + \Lambda) \Delta_n + \Phi^{1/2} \Delta B_n,
#' }{
#' Y_n+1 = Y_n + (\Gamma Y_n + \Lambda) \Delta_n + \Phi^(1/2) \Delta B_n,
#' }
#' where \eqn{Y_n = Y(t_n)}, \eqn{\Delta_n = t_{n+1} - t_n}{\Delta_n = t_n+1 - t_n} and
#' \deqn{
#' \Delta B_n = B(t_{n+1}) - B(t_n) \stackrel{\textnormal{ind}}{\sim} \mathcal N(0, \Delta_n).
#' }{
#' \Delta B_n = B(t_n+1) - B(t_n) \sim iid N(0, \Delta_n).
#' }
#' Thus, \eqn{Y_0, \ldots, Y_N} is multivariate normal Markov chain for which the marginal distribution of any subset of timepoints and/or components can be efficiently calculated using the Kalman filter.  This can be used to check the MCMC output of [sde.post()] as in the example.
#'
#' @example examples/mou.loglik.R
#' @export
mou.loglik <- function(X, dt, nvar.obs, Gamma, Lambda, Phi, mu0, Sigma0) {
  nobs <- nrow(X)
  ndims <- ncol(X)
  # check parameters
  if(!identical(dim(Gamma), c(ndims, ndims))) {
    stop("X and Gamma have incompatible dimensions.")
  }
  if(!is.vector(Lambda) || length(Lambda) != ndims) {
    stop("X and Lambda have incompatible dimensions.")
  }
  if(!identical(dim(Phi), c(ndims, ndims))) {
    stop("X and Phi have incompatible dimensions.")
  }
  if(any(diag(chol(Phi)) <= 0)) {
    stop("Phi must be a variance matrix.")
  }
  # check dt and nvar.obs
  if(length(dt) == 1) dt <- rep(dt, nobs - 1)
  if(length(dt) != nobs - 1) stop("x and dt have incompatible sizes.")
  if(missing(nvar.obs)) nvar.obs <- ndims
  if(length(nvar.obs) == 1) nvar.obs <- rep(nvar.obs, nobs)
  if(length(nvar.obs) != nobs) stop("nvar.obs and x have incompatible sizes.")
  if(!all(nvar.obs == round(nvar.obs)) || any(nvar.obs < 0) ||
     any(nvar.obs > ndims)) {
    stop("nvar.obs must be an integer scalar or vector with values between 0 and n.")
  }
  lpdf <- rep(NA, nobs)
  if(missing(mu0)) mu0 <- rep(0, ndims)
  if(missing(Sigma0)) Sigma0 <- diag(ndims)
  # initialize
  mu_n <- mu0
  Sigma_n <- Sigma0
  qn <- nvar.obs[1] # evaluate PDF
  lpdf[1] <- .lmvn(X[1,1:qn] - mu_n[1:qn], Sigma_n[1:qn,1:qn,drop=FALSE])
  if(nobs == 1) return(lpdf)
  # recursions
  for(nn in 1:(nobs-1)) {
    # update
    A_n <- diag(ndims) + Gamma * dt[nn]
    b_n <- Lambda * dt[nn]
    C_n <- Phi * dt[nn]
    if(qn == ndims) {
      # all elements of Y_n are observed
      mu_n <- A_n %*% X[nn,] + b_n
      Sigma_n <- C_n
    } else if(qn == 0) {
      # no elements of Y_n are observed
      mu_n <- A_n %*% mu_n + b_n
      Sigma_n <- C_n + A_n %*% Sigma_n %*% t(A_n)
    } else {
      # some but not all elements of Y_n are observed
      condMV <- .cmvn(X[nn,], mu_n, Sigma_n, qn)
      m_n <- condMV$m
      V_n <- condMV$V
      mu_n <- b_n + A_n %*% c(X[nn,1:qn],m_n)
      Sigma_n <- C_n + A_n[,(qn+1):ndims,drop=FALSE] %*%
        V_n %*% t(A_n[,(qn+1):ndims,drop=FALSE])
    }
    # evaluate new PDF
    qn <- nvar.obs[nn+1]
    lpdf[nn+1] <- .lmvn(X[nn+1,1:qn] - mu_n[1:qn],
                        Sigma_n[1:qn,1:qn,drop=FALSE])
  }
  sum(lpdf)
}

# helper functions

# solve method for variance matrices
# optionally computes log determinant as well
.solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}

# log-pdf calculation for N(0, Sigma), one observation only
.lmvn <- function(x, Sigma) {
  if(length(x) == 0) return(0)
  l2pi <- 1.837877066409345483560659472811 # log(2*pi)
  iV <- .solveV(Sigma, x, ldV = TRUE)
  -.5 * (crossprod(x, iV$y)[1] + iV$ldV + length(x) * l2pi)
}

# conditional mean and variance of a multivariate normal
# x, mu, Sigma are fully observed.
# qn is between 1 and length(x)
.cmvn <- function(x, mu, Sigma, qn) {
  ndims <- length(x)
  rn <- ndims-qn
  x <- x[1:qn] # selecting the observed Xn's
  m2 <- mu[1:qn]
  m1 <- mu[(qn+1):ndims]
  # double check the elements selection
  V22 <- Sigma[1:qn,1:qn,drop=FALSE] # qn x qn
  V21 <- Sigma[1:qn,(qn+1):ndims,drop=FALSE] # qn x (d - qn)
  V12 <- Sigma[(qn+1):ndims,1:qn,drop=FALSE] # (d - qn) x qn
  V11 <- Sigma[(qn+1):ndims,(qn+1):ndims,drop=FALSE] # (d - qn) x (d - qn)
  # all inversion in one step
  iV <- crossprod(V21, .solveV(V22, cbind(V21, x - m2)))
  mn <- m1 + iV[,rn+1]
  Vn <- V11 - iV[,1:rn]
  list(m = mn, V = Vn)
}
