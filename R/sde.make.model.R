#' Create an SDE model object.
#'
#' Compiles the C++ code for various SDE-related algorithms and makes the routines available within R.
#'
#' @param ModelFile Path to the header file where the SDE model is defined.
#' @param PriorFile Path to the header file where the SDE prior is defined.  See [sde.prior()] for details.
#' @param data.names Vector of names for the SDE components.  Defaults to `X1,...,Xd`.
#' @param param.names Vector of names for the SDE parameters.  Defaults to `theta1,...,thetap`.
#' @param hyper.check A function with arguments `hyper`, `param.names`, and `data.names` used for passing the model hyper parameters to the C++ code.  See [mvn.hyper.check()] for details.
#' @param OpenMP Logical; whether the model is compiled with `OpenMP` for C++ level parallelization.
#' @param ... additional arguments to [Rcpp::sourceCpp()] for compiling the C++ code.
#'@return An `sde.model` object, consisting of a list with the following elements:
#' \describe{
#' \item{`ptr`}{Pointer to C++ sde object (`sdeRobj`) implementing the member functions: drift/diffusion, data/parameter validators, loglikelihood, prior distribution, forward simulation, MCMC algorithm for Bayesian inference.}
#' \item{`ndims`, `nparams`}{The number of SDE components and parameters.}
#' \item{`data.names`, `param.names`}{The names of the SDE components and parameters.}
#' \item{`omp`}{A logical flag for whether or not the model was compiled for multicore functionality with `OpenMP`.}
#' }
#' @seealso [sde.drift()], [sde.diff()], [sde.valid()], [sde.loglik()], [sde.prior()], [sde.sim()], [sde.post()].
#' @example examples/sde.make.model.R
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
  ## # save all sdeObj pointers in the package environment
  ## # as sdeObj pointers don't gc properly when R object is overwritten
  ## globalptr <- gsub("^.", "", tempfile(pattern = "sdeObj_", tmpdir = ""))
  ## # compile C++ code
  ## cppFile <- .copy.cpp.files(ModelFile, PriorFile, OpenMP)
  ## if(OpenMP) old.env <- .omp.set()
  ## sourceCpp(cppFile, env = .msdeglobalenv, ...)
  ## if(OpenMP) .omp.unset(env = old.env)
  ## assign(globalptr, .msdeglobalenv$.sde_MakeModel(), envir = .msdeglobalenv)
  ## sptr <- .msdeglobalenv[[globalptr]]
  # compile C++ code
  cppFile <- .copy.cpp.files(ModelFile, PriorFile, OpenMP)
  if(OpenMP) old.env <- .omp.set()
  sourceCpp(cppFile, ...)
  if(OpenMP) .omp.unset(env = old.env)
  cobj <- eval(parse(
    text = paste0("new(",
                  gsub("_Module[.]cpp$", "", basename(cppFile)),
                  ")")
  ))
  # extract ndims and nparams
  ndims <- cobj$nDims()
  nparams <- cobj$nParams()
  # parameter and data names
  if(missing(data.names)) data.names <- paste0("X", 1:ndims)
  if(missing(param.names)) param.names <- paste0("theta", 1:nparams)
  if(length(data.names) != ndims) {
    stop("data.names has wrong length.")
  }
  if(length(param.names) != nparams) {
    stop("param.names has wrong length.")
  }
  sde.model <- list(cobj = cobj, ndims = ndims, nparams = nparams,
                    data.names = data.names, param.names = param.names,
                    hyper.check = hyper.check, omp = OpenMP)
  # output
  class(sde.model) <- "sde.model"
  sde.model
}

#--- helper functions for keeping track of models ------------------------------

#' Global environment for tracking model compilation.
#'
#' @details This environment contains a global variable called `models` consisting of a `data.frame` where each row is a model and the columns are:
#' \describe{
#'   \item{`id`}{A unique model identifier.}
#'   \item{`sde`}{The md5sum of the given `ModelFile`.}
#'   \item{`prior`}{The md5sum of the given `PriorFile`.}
#'   \item{`omp`}{The logical flag given by `OpenMP`.}
#' }
#' This is needed because `sde.make.model()` returns an instance of an object rather than a class definition that can be used to instantiate objects at will.  So, if `sde.make.model()` is given identical C++ files `ModelFile` and `PriorFile`, then we don't want to recompile the C++ code but rather use the existing class definition to instantiate the return object.
#' @noRd
.msdeglobalenv <- new.env(parent = globalenv())

#' Obtain the model ID from C++ code specification.
#'
#' @param ModelFile Path to the header file where the SDE model is defined.
#' @param PriorFile Path to the header file where the SDE prior is defined.
#' @param OpenMP Logical; whether the model is compiled with `OpenMP` for C++ level parallelization.
#'
#' @return A unique ID corresponding to the model.  The ID is added to the list stored in `.msdeglobalenv$models` which contains all unique models so far, i.e., for which `md5sum(ModelFile)`, `md5sum(PriorFile)`, and `OpenMP` are unique.
#'
#' @noRd
.cpp.model.id <- function(ModelFile, PriorFile, OpenMP) {
  mod <- data.frame(id = tempfile(pattern = "msde_"),
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

#' Create a tmp directory containing model files.
#'
#' @param ModelFile Path to the header file where the SDE model is defined.
#' @param PriorFile Path to the header file where the SDE prior is defined.
#' @param OpenMP Logical; whether the model is compiled with `OpenMP` for C++ level parallelization.
#'
#' @return The full path to the `cpp` file to be compiled, i.e.,  inside `tempdir()`.  As a side effect, copies `ModelFile` and `PriorFile` to `tempdir()`.
#'
#' @details Compilation of identical files is avoided through \pkg{Rcpp}'s caching mechanism which doesn't recompile a `cpp` file with the same name and md5sum.  So, if `.cpp.model.id()` determines that the input files and the `OpenMP` flag are identical, it will give the `cpp` file to be compiled the existing name from `.msdeglobalenv$models`, and \pkg{Rcpp} will do the rest.
#' @noRd
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
  fname <- .cpp.model.id(ModelFile, PriorFile, OpenMP)
  model.name <- basename(fname)
  fname <- paste0(fname, "_Module.cpp")
  # create file from template
  # FIXME: can potentially do this with Rcpp::exposeClass...
  template <- readLines(file.path(.msde_include_path,
                                  "templates", "sdeModule.cpp"))
  template.data <- list(
    ExportFile = paste0(model.name, "_Module.cpp"),
    ModelFile = "sdeModel.h",
    PriorFile = "sdePrior.h",
    sdeModel = "sdeModel",
    sdePrior = "sdePrior",
    ModuleName = model.name,
    RClassName = model.name
  )
  module <- whisker.render(template = template, data = template.data)
  cat(module, file = fname, sep = "\n")
  ## fname <- paste0(.cpp.model.id(ModelFile, PriorFile, OpenMP), "_Exports.cpp")
  ## file.copy(from = file.path(.msde_tools_path, "sdeMakeModel.cpp"),
  ##           to = fname,
  ##           overwrite = TRUE, copy.date = TRUE)
  fname
}

#' Expose the `sdeRobj<sdeModel, sdePrior>` class to R.
#'
#' @param ModuleFile The name of the file (or connection) to which to write the \pkg{Rcpp} module.  See `Rcpp::exposeClass()`.
#' @param ModuleName The name of the module containing the class.
#' @param RClassName The name of the class on the R side.
#' @param sdeModel,sdePrior The C++ class names for the SDE model and prior.
#' @param ModelIncludes A vector of strings, each of which is a line of C++ code, providing what's necessary to include the definitions of `sdeModel` and `sdePrior`.  The default value is `c('#include "sdeModel.h"', '#include "sdePrior.h"')`.
#' @param RFile Optional R file to wrap the RC class on the R side.  See `Rcpp::exposeClass()`.
#'
#' @return Nothing.  Called for the side effect of creating `ModuleFile` and (optionally) `RFile`.
#'
#' @noRd
.sde.expose.model <- function(ModuleFile, ModuleName, RClassName,
                              sdeModel, sdePrior,
                              ModelIncludes, RFile = FALSE) {
  CppClassName <- paste0("sdeRobj<", sdeModel, ", ", sdePrior, "> ")
  ClassTypedef <- gsub("[^a-zA-Z0-9]+", "_",
                       paste0('sdeRobj_', sdeModel, '_', sdePrior))
  # header includes
  msdeIncludes <- c('//[[Rcpp::depends(msde)]]',
                    '//[[Rcpp::depends(RcppProgress)]]',
                    '#include <sdeRobj.h>')
  if(missing(ModelIncludes)) {
    ModelIncludes <- c('#include "sdeModel.h"',
                       '#include "sdePrior.h"')

  }
  HeaderIncludes <- paste(c(
    msdeIncludes,
    ModelIncludes,
    paste0('typedef ', CppClassName, ' ', ClassTypedef, ';')
  ), sep = "\n")
  Rcpp::exposeClass(class = RClassName,
                    constructors = list(""),
                    methods = c("get_nDims", "get_nParams",
                                "isData", "isParams", "Drift", "Diff",
                                "LogLik", "Prior", "Sim", "Post"),
                    header = HeaderIncludes,
                    fields = character(),
                    ## CppClass = ClassTypedef,
                    CppClass = CppClassName,
                    module = ModuleName,
                    rename = c(nDims = "get_nDims", nParams = "get_nParams",
                               Loglik = "LogLik"),
                    file = ModuleFile, Rfile = RFile)
}

#--- omp set and unset ---------------------------------------------------------

#' Adds `-fopenmp` flags to `PKG_CXXFLAGS` and `PKG_LIBS`.
#'
#' @noRd
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

#' Resets environment variables modified by `.omp.set()` to those supplied in `env`.
#'
#' @noRd
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

