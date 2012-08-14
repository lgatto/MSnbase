.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
                        paste("\nThis is MSnbase version",packageVersion("MSnbase"),"\n",
                              " Read '?MSnbase' and references therein for information\n",
                              " about the package and how getting started.\n"))
  addVigs2WinMenu("MSnbase-demo")
  addVigs2WinMenu("MSnbase-dev")
}

