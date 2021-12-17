.onLoad <- function(libname, pkgname) {
  # path to msde header-only library
  assign(".msde_include_path",
         file.path(libname, pkgname, "include"),
         envir = parent.env(environment()))
  # path to R-specific on-the-fly compile mechanism
  assign(".msde_tools_path",
         file.path(libname, pkgname, "include", "R"),
         envir = parent.env(environment()))
  # Pre-compiled sde models
  loadModule("class_msde_hestModel", TRUE)
  loadModule("class_msde_lotvolModel", TRUE)
  loadModule("class_msde_biouModel", TRUE)
  loadModule("class_msde_pgnetModel", TRUE)
  loadModule("class_msde_eouModel", TRUE)
}
