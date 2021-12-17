#' R6 class wrapping `sdeRobj<{{sdeModel}}, {{sdePrior}}>`.
#'
#' @details Should this class be hidden?
#'
#' If yes, the main advantage is that the class does not need to be documented.  It is then the instantiated object which would potentially be exported, since it would be called with the `sde.{method}` generics provided by **msde**.
#'
#' However, if one wrote a package that provided SDE models, the package would probably want to export the class constructor, so as to instantiate multiple model objects.
#'
#' So let's do the following for now:
#'
#' Let's create a new separate class for each polymorphic object, e.g., `msde_MySdeModel`.  We'll put all the class definitions into an R file.
#'
#' The class definitions are not exported from the package.  In this particular case, they are not meant to be used directly, always from `sde.{method}`.  However, this doesn't resolve the issue of exporting models from users' packages.  I think it's easiest for now to not provide this feature.
#' @noRd
{{RClassName}} <- R6::R6Class(

  classname = "{{RClassName}}",

  private = list(
    .ptr = NULL,
    .ndims = NULL,
    .nparams = NULL,
    .data.names = NULL,
    .param.names = NULL,
    .omp = NULL,
    .isData = {{RClassName}}_isData,
    .isParams = {{RClassName}}_isParams,
    .Drift = {{RClassName}}_Drift,
    .Diff = {{RClassName}}_Diff,
    .Loglik = {{RClassName}}_Loglik,
    .Prior = {{RClassName}}_Prior,
    .Sim = {{RClassName}}_Sim,
    .Post = {{RClassName}}_Post
  ),

  active = list(

    #' @field ptr Pointer to C++ `sdeRobj<{{sdeModel}}, {{sdePrior}}>`.
    ptr = function(value) {
      if(missing(value)) {
        private$.ptr
      } else {
        stop("$ptr is read-only.", call. = FALSE)
      }
    },

    #' @field ndim Number of SDE components.
    ndims = function(value) {
      if(missing(value)) {
        private$.ndims
      } else {
        stop("$ndims is read-only.", call. = FALSE)
      }
    },

    #' @field nparams Number of model parameters.
    nparams = function(value) {
      if(missing(value)) {
        private$.nparams
      } else {
        stop("$nparams is read-only.", call. = FALSE)
      }
    },

    #' @field data.names Names of the SDE components.
    data.names = function(value) {
      if(missing(value)) {
        private$.data.names
      } else {
        if(length(value) != self$ndims) {
          stop("$data.names has wrong length.")
        } else {
          private$.data.names <- value
        }
      }
    },

    #' @field param.names Names of the model parameters.
    param.names = function(value) {
      if(missing(value)) {
        private$.param.names
      } else {
        if(length(value) != self$nparams) {
          stop("$param.names has wrong length.")
        } else {
          private$.param.names <- value
        }
      }
    },

    #' @field omp A logical flag for whether or not the model was compiled for multicore functionality with \code{OpenMP}.
    omp = function(value) {
      if(missing(value)) {
        private$.omp
      } else {
        stop("$omp is read-only.", call. = FALSE)
      }
    }
  ),

  public = list(
    isData = function(...) private$.isData(private$.ptr, ...),
    isParams = function(...) private$.isParams(private$.ptr, ...),
    Drift = function(...) private$.Drift(private$.ptr, ...),
    Diff = function(...) private$.Diff(private$.ptr, ...),
    Loglik = function(...) private$.Loglik(private$.ptr, ...),
    Prior = function(...) private$.Prior(private$.ptr, ...),
    Sim = function(...) private$.Sim(private$.ptr, ...),
    Post = function(...) private$.Diff(private$.ptr, ...),
    hyper.check = NULL,

    initialize = function(data.names, param.names,
                          hyper.check,
                          OpenMP = FALSE) {
      # initialize sdeRobj<{{sdeModel}}, {{sdePrior}}> object in C++
      private$.ptr <- {{RClassName}}_ctor()
      # set ndims and nparams
      private$.ndims <- {{RClassName}}_nDims(private$.ptr)
      private$.nparams <- {{RClassName}}_nParams(private$.ptr)
      # parameter and data names
      if(missing(data.names)) data.names <- paste0("X", 1:self$ndims)
      if(missing(param.names)) param.names <- paste0("theta", 1:self$nparams)
      self$data.names <- data.names
      self$param.names <- param.names
      # remaining members
      self$hyper.check <- hyper.check
      private$.omp <- OpenMP
    }
  )
)
