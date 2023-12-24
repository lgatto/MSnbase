.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
      paste("\nThis is MSnbase version", packageVersion("MSnbase"), "\n",
            " Visit https://lgatto.github.io/MSnbase/ to get started.\n",
            "Consider switching to the 'R for Mass Spectrometry'\n",
            "packages - see https://RforMassSpectrometry.org for details.\n"))
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
    fast_load <- TRUE
    ## Disable "fast reading" on macOS. See issue #170
    ## Enabling fast reading also on macOS to avoid spurious errors.
    ## if (Sys.info()["sysname"] == "Darwin")
    ##     fast_load <- FALSE
    msOps <- list(PARALLEL_THRESH = 1000,
                  fastLoad = fast_load,      # disable reading header before peaks.
                  ## sortMethod = sortMeth,
                  verbose = FALSE)
    options(MSnbase = msOps)
}
