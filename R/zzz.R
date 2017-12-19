.onLoad <- function(libname, pkgname) {
  # path to msde header-only library
  assign(".msde_include_path",
         file.path(libname, pkgname, "include"),
         envir = parent.env(environment()))
  # path to R-specific on-the-fly compile mechanism
  assign(".msde_tools_path",
         file.path(libname, pkgname, "include", "R"),
         envir = parent.env(environment()))
}
