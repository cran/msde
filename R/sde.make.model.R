#' Create an SDE model object.
#'
#' Compiles the C++ code for various SDE-related algorithms and makes the routines available within R.
#' @param ModelFile Path to the header file where the SDE model is defined.
#' @param PriorFile Path to the header file where the SDE prior is defined.  See \code{\link{sde.prior}} for details.
#' @param data.names Vector of names for the SDE components.  Defaults to \code{X1,...,Xd}.
#' @param param.names Vector of names for the SDE parameters.  Defaults to \code{theta1,...,thetap}.
#' @param hyper.check A function with arguments \code{hyper}, \code{param.names}, and \code{data.names} used for passing the model hyper parameters to the C++ code.  See \code{\link{mvn.hyper.check}} for details.
#' @param OpenMP Logical; whether the model is compiled with \code{OpenMP} for C++ level parallelization.
#' @param ... additional arguments to \code{Rcpp::sourceCpp} for compiling the C++ code.
#'@return An \code{sde.model} object, consisting of a list with the following elements:
#' \describe{
#' \item{\code{model.ptr}}{Pointer to C++ sde object (\code{sdeCobj}) implementing the member functions: drift/diffusion, data/parameter validators, loglikelihood, prior distribution, forward simulation, MCMC algorithm for Bayesian inference.}
#' \item{\code{ndims, nparams}}{The number of SDE components and parameters.}
#' \item{\code{data.names, param.names}}{The names of the SDE components and parameters.}
#' \item{\code{omp}}{A logical flag for whether or not the model was compiled for multicore functionality with \code{OpenMP}.}
#' }
#' @seealso \code{\link{sde.drift}}, \code{\link{sde.diff}}, \code{\link{sde.valid}}, \code{\link{sde.loglik}}, \code{\link{sde.prior}}, \code{\link{sde.sim}}, \code{\link{sde.post}}.
#' @importFrom Rcpp sourceCpp
#' @importFrom methods formalArgs
#' @importFrom tools md5sum
#' @examples
#' # header (C++) file for Heston's model
#' hfile <- sde.examples("hest", file.only = TRUE)
#' cat(readLines(hfile), sep = "\n")
#'
#' \donttest{
#' # compile the model
#' param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
#' data.names <- c("X", "Z")
#' hmod <- sde.make.model(ModelFile = hfile,
#'                        param.names = param.names,
#'                        data.names = data.names)
#'
#' hmod
#' }
#' @export
sde.make.model <- function(ModelFile, PriorFile = "default",
                           data.names, param.names, hyper.check,
                           OpenMP = FALSE, ...) {
  # prior specification
  if(PriorFile == "default") {
    PriorFile <- file.path(.msde_include_path, "mvnPrior.h")
    if(!missing(hyper.check)) {
      warning("Custom hyper.check ignored for default prior.")
    }
    hyper.check <- mvn.hyper.check
  } else {
    if(missing(hyper.check)) {
      stop("Must provide hyper.check for custom prior.")
    }
    if(!identical(formalArgs(hyper.check),
                  c("hyper", "param.names", "data.names"))) {
      stop("hyper.check must have formal arguments: hyper, param.names, data.names.")
    }
  }
  # save all sdeObj pointers in the package environment
  # as sdeObj pointers don't gc properly when R object is overwritten
  globalptr <- gsub("^.", "", tempfile(pattern = "sdeObj_", tmpdir = ""))
  # compile C++ code
  cppFile <- .copy.cpp.files(ModelFile, PriorFile, OpenMP)
  if(OpenMP) old.env <- .omp.set()
  sourceCpp(cppFile, env = .msdeglobalenv, ...)
  if(OpenMP) .omp.unset(env = old.env)
  assign(globalptr, .msdeglobalenv$.sde_MakeModel(), envir = .msdeglobalenv)
  sptr <- .msdeglobalenv[[globalptr]]
  # extract ndims and nparams
  ndims <- .sde_nDims(sptr)
  nparams <- .sde_nParams(sptr)
  # parameter and data names
  if(missing(data.names)) data.names <- paste0("X", 1:ndims)
  if(missing(param.names)) param.names <- paste0("theta", 1:nparams)
  if(length(data.names) != ndims) {
    stop("data.names has wrong length.")
  }
  if(length(param.names) != nparams) {
    stop("param.names has wrong length.")
  }
  sde.model <- list(ptr = sptr, ndims = ndims, nparams = nparams,
                    data.names = data.names, param.names = param.names,
                    hyper.check = hyper.check, omp = OpenMP)
  # output
  class(sde.model) <- "sde.model"
  sde.model
}

#--- keep track of original models ---------------------------------------------

# global variable: md5sum of model/prior pairs, modelID
.msdeglobalenv <- new.env(parent = globalenv())

# outputs the file "id"
.cpp.model.id <- function(ModelFile, PriorFile, OpenMP) {
  mod <- data.frame(id = tempfile(pattern = "msde-"),
                    sde = md5sum(ModelFile)[1],
                    prior = md5sum(PriorFile)[1],
                    omp = OpenMP,
                    stringsAsFactors = FALSE)
  models <- .msdeglobalenv$models
  if(is.null(models)) models <- mod
  same <- (mod$sde == models$sde)
  same <- same & (mod$prior == models$prior)
  same <- same & (mod$omp == models$omp)
  if(any(same)) {
    mod$id <- models$id[which(same)[1]]
  } else {
    models <- rbind(models, mod)
  }
  # save environment variable
  assign(x = "models", value = models, envir = .msdeglobalenv)
  mod$id
}

.copy.cpp.files <- function(ModelFile, PriorFile, OpenMP) {
  # prior file
  fname <- file.path(tempdir(), "sdePrior.h")
  flag <- file.copy(from = PriorFile,
                    to = fname,
                    overwrite = TRUE, copy.date = TRUE)
  if(!flag) {
    stop("PriorFile \"", PriorFile, "\" not found.")
  }
  # model file
  fname <- file.path(tempdir(), "sdeModel.h")
  flag <- file.copy(from = ModelFile,
                    to = fname,
                    overwrite = TRUE, copy.date = TRUE)
  if(!flag) {
    stop("ModelFile \"", ModelFile, "\" not found.")
  }
  # export file
  fname <- paste0(.cpp.model.id(ModelFile, PriorFile, OpenMP), "_Exports.cpp")
  file.copy(from = file.path(.msde_tools_path, "sdeMakeModel.cpp"),
            to = fname,
            overwrite = TRUE, copy.date = TRUE)
  fname
}

#--- omp set and unset ---------------------------------------------------------

# adds -fopenmp flags to PKG_CXXFLAGS and PKG_LIBS

.omp.set <- function() {
  cxx <- Sys.getenv(x = "PKG_CXXFLAGS", unset = NA)
  libs <- Sys.getenv(x = "PKG_LIBS", unset = NA)
  env <- list(cxx = cxx, libs = libs)
  Sys.setenv(PKG_CXXFLAGS = ifelse(is.na(cxx),
                                "-fopenmp", paste("-fopenmp", cxx)))
  Sys.setenv(PKG_LIBS = ifelse(is.na(libs),
                               "-fopenmp", paste("-fopenmp", libs)))
  env
}

.omp.unset <- function(env) {
  if(is.na(env$cxx)) {
    Sys.unsetenv(x = "PKG_CXXFLAGS")
  } else {
    Sys.setenv(PKG_CXXFLAGS = env$cxx)
  }
  if(is.na(env$libs)) {
    Sys.unsetenv(x = "PKG_LIBS")
  } else {
    Sys.setenv(PKG_LIBS = env$libs)
  }
}

