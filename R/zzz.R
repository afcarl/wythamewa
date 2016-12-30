.onLoad <- function(libname, pkgname) { }

.onAttach <- function(...) {
  theLib <- dirname(system.file(package = "wythamewa"))
  pkgdesc <- packageDescription("wythamewa", lib.loc = theLib)
  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
  msg <- paste("wythamewa (Version ", pkgdesc$Version, ")", sep = "")
  msg <- paste(msg,"Type ?wythamewa for a quick introduction.",sep="\n")
  packageStartupMessage(msg)
}