.onLoad <- function(libname, pkgname) {
  assign(".msde_tools_path",
         file.path(libname, pkgname, "tools"),
         envir = parent.env(environment()))
  assign(".msde_include_path",
         file.path(libname, pkgname, "include"),
         envir = parent.env(environment()))
}
