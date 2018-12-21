.onLoad <- function(libname, pkgname) {
  rJava::.jpackage(pkgname, lib.loc = libname)
}

.onAttach <- function(libname, pkgname){
  msg <- paste0(
    "rTASSEL (experimental)\n",
    "Use with caution!\n"
  )
  packageStartupMessage(msg)
}