.onLoad <- function(libname, pkgname) {
  # Package loading function
  invisible()
}

.onUnload <- function(libpath) {
  library.dynam.unload("TWoLife", libpath)
}