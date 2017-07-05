.onLoad <- function(libname, pkgname) {
  assign(".msde_src_path", file.path(libname, pkgname, "msde", "src"),
         envir = parent.env(environment()))
  assign(".msde_examples_path", file.path(libname, pkgname, "msde", "examples"),
         envir = parent.env(environment()))
}
