#' Example SDE models.
#'
#' Provides sample `C++` code for several SDE models.
#'
#' @param model Character string giving the name of a sample model.  Possible values are: `hest`, `pgnet`, `lotvol`, `biou`, `eou`.  See Details.
#' @param file.only If `TRUE` returns only the path to the header file containing the `sdeModel` object implementation.
#'
#' @return An `sde.model` object, or the path to the C++ model header file.
#'
#' @details All pre-compiled models are with the default prior and with `OpenMP` disabled.  A full description of the example models can be found in the package vignette; to view it run `vignette("msde-exmodels")`.
#'
#' @seealso [sde.make.model()] for `sde.model` objects, [mvn.hyper.check()] for specification of the default prior.
#' @example examples/sde.examples.R
#' @export
sde.examples <- function(model = c("hest", "pgnet", "lotvol", "biou", "eou"),
                         file.only = FALSE) {
  model <- match.arg(model)
  if(model == "hest") {
    ModelFile <- file.path(.msde_include_path, "hestModel.h")
    param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
    data.names <- c("X", "Z")
    cobj <- new(msde_hestModel)
    ## cobj <- msde_hestModel()
  } else if(model == "pgnet") {
    ModelFile <- file.path(.msde_include_path, "pgnetModel.h")
    param.names <- paste0("theta", 1:8)
    data.names <- c("R", "P", "Q", "D")
    cobj <- new(msde_pgnetModel)
    ## cobj <- msde_pgnetModel()
  } else if(model == "biou") {
    ModelFile <- file.path(.msde_include_path, "biouModel.h")
    param.names <- c("Gamma11", "Gamma21", "Gamma12", "Gamma22",
                     "Lambda1", "Lambda2",
                     "Psi11", "Psi12", "Psi22")
    data.names <- c("Y1","Y2")
    cobj <- new(msde_biouModel)
    ## cobj <- msde_biouModel()
  } else if(model == "lotvol") {
    ModelFile <- file.path(.msde_include_path, "lotvolModel.h")
    param.names <- c("alpha", "beta", "gamma")
    data.names <- c("H", "L")
    cobj <- new(msde_lotvolModel)
    ## cobj <- msde_lotvolModel()
  } else if(model == "eou") {
    ModelFile <- file.path(.msde_include_path, "eouModel.h")
    param.names <- c("alpha", "gamma", "eta", "sigma", "rho")
    data.names <- c("X", "V")
    cobj <- new(msde_eouModel)
    ## cobj <- msde_eouModel()
  }
  if(!file.only) {
    sde.model <- list(cobj = cobj,
                      ndims = length(data.names), nparams = length(param.names),
                      data.names = data.names, param.names = param.names,
                      hyper.check = mvn.hyper.check, omp = FALSE)
    class(sde.model) <- "sde.model"
  } else sde.model <- ModelFile
  sde.model
}
