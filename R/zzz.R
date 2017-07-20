.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
      paste("\nThis is MSnbase version", packageVersion("MSnbase"), "\n",
            " Read '?MSnbase' and references therein for information\n",
            " about the package and how to get started.\n"))

  if (interactive() && .Platform$OS.type == "windows" &&
      .Platform$GUI == "Rgui") {
      Biobase::addVigs2WinMenu("MSnbase")
  }
}

.onLoad <- function(libname, pkgname) {
    ## Add MSnbase options.
    ## Use radix sorting for R >= 3.3
    ## sortMeth <- "auto"
    ## if (as.numeric(R.Version()$major) >= 3 & as.numeric(R.Version()$minor) >= 3)
    ##     sortMeth <- "radix"
    msOps <- list(PARALLEL_THRESH = 1000,
                  fastLoad = TRUE,      # disable reading header before peaks.
                  ## sortMethod = sortMeth,
                  verbose = FALSE)
    options(MSnbase = msOps)
}
