#--- argument checking ---------------------------------------------------------

# check if vars are in var.names, non-NULL, no duplicates
is.valid.vars <- function(vars, var.names) {
  valid <- !is.null(vars)
  valid <- valid && !anyDuplicated(vars)
  valid <- valid && all(vars %in% var.names)
  valid
}

# check if hyperparameters are a list of NULL or vector-doubles
is.valid.hyper <- function(phi) {
  valid <- is.list(phi)
  valid <- valid && all(sapply(phi, function(x) {
    is.null(x) || (is.double(x) && is.vector(x))
  }))
  valid
}

# check a vector of lengths and make sure that they
# are all the same or equal to 1
is.valid.nreps <- function(nreps) {
  nreps <- unique(nreps)
  nreps <- c(1, nreps[nreps > 1])
  length(nreps) <= 2
}

# format data for sde.drift, sde.diff, sde.loglik
# output is a matrix or 3-d array with dimensions 1 and 2 permuted
.format.data <- function(x, data.names, type = c("matrix", "array"),
                         strict = FALSE, debug = FALSE) {
  ndims <- length(data.names)
  type <- match.arg(type)
  if(debug) browser()
  # check is.numeric, ndims, and data.names
  if(!is.numeric(x) || !all(is.finite(x))) {
    stop("x must be numeric with finite values.")
  }
  if(type == "matrix") {
    if(is.vector(x)) {
      if(strict) {
        stop("x must be a matrix.")
      } else {
        x <- t(x)
      }
    }
    if(!is.matrix(x)) {
      stop("x must be a vector or matrix.")
    }
    if(ncol(x) != ndims) {
      stop("dimensions of x are incompatible with ndims.")
    }
    if(!is.null(colnames(x)) && !identical(colnames(x), data.names)) {
      stop("names of x do not match data.names.")
    }
    x <- t(x)
  } else {
    if(is.matrix(x)) {
      if(strict) {
        stop("x must be an array.")
      } else {
        dn <- dimnames(x)
        x <- array(x, dim = c(dim(x), 1))
        if(!is.null(dn)) dimnames(x) <- c(dn, list(NULL))
      }
    }
    if(length(dim(x)) != 3) {
      stop("x must be a matrix or 3-d array.")
    }
    if(dim(x)[2] != ndims) {
      stop("dimensions of x are incompatible with ndims.")
    }
    if(!is.null(dimnames(x)) && !identical(dimnames(x)[[2]], data.names)) {
      stop("names of x do not match data.names.")
    }
    x <- aperm(x, c(2,1,3))

  }
  if(is.null(dimnames(x))) dimnames(x) <- rep(list(NULL), length(dim(x)))
  dimnames(x)[[1]] <- data.names
  x
}

# format parameters for sde.drift, sde.diff, sde.loglik
# output is a matrix with dimensions 1 and 2 permuted
.format.params <- function(theta, param.names) {
  nparams <- length(param.names)
  if(is.vector(theta)) theta <- t(theta)
  # check is.numeric, nparams, param.names
  if(!is.numeric(theta) || !all(is.finite(theta))) {
    stop("theta must be a numeric with finite values.")
  }
  if(!is.matrix(theta)) {
    stop("theta must be a vector or matrix.")
  }
  if(ncol(theta) != nparams) {
    stop("dimensions of theta are incompatible with nparams.")
  }
  if(!is.null(colnames(theta)) && !identical(colnames(theta), param.names)) {
    stop("names of theta do not match param.names.")
  }
  if(is.null(dimnames(theta))) {
    dimnames(theta) <- rep(list(NULL), length(dim(theta)))
  }
  dimnames(theta)[[2]] <- param.names
  t(theta)
}

## # check that variable names are compatible with param.names and data.names
## # var.names adds information to error message since .check.vars is not
## # an exported function
## .check.vars <- function(vars, param.names, data.names,
##                         vnames = "vars") {
##   nparams <- length(param.names)
##   ndims <- length(data.names)
##   var.names <- c(param.names, data.names)
##   if(is.null(vars) || !is.character(vars)) {
##     stop(vnames, " must be a non-null character vector.")
##   }
##   if(anyDuplicated(vars)) {
##     stop(vnames, " must be unique.")
##   }
##   if(!all(vars %in% var.names)) {
##     stop(vnames, " not in param.names or data.names.")
##   }
##   ## vars2 <- c(param.names[!fixed.params], data.names[ndims:1 <= nmiss])
##   ## var.id <- var.names %in% vars
##   ## if(!all(var.id) && !all(var.id == (var.names %in% vars2))) {
##   ##   stop("vars and (fixed.params,nmiss) specify different active sets.")
##   ## }
##   TRUE
## }

## # check that phi is a list of double vectors
## .check.hyper <- function(phi, pname = "phi") {
##   if(!is.list(phi)) stop(pname, " must be a list.")
##   is.valid.phi <- sapply(phi, function(x) {
##     is.null(x) || (is.double(x) & is.vector(x))
##   })
##   is.valid.phi <- all(is.valid.phi)
##   if(!is.valid.phi) {
##     stop("Each element of ", pname, " must be NULL or a double vector.")
##   }
##   TRUE
## }
