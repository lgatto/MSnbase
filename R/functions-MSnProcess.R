show.MSnProcess <- function(object) {
  cat("- - - Processing information - - -\n")
  for (proc in object@processing)
    cat(" ",proc,"\n")
  cat(" MSnbase version:",packageVersion("MSnbase"),"\n")
  cat(" Xcms version:",object@xcmsVersion,"\n")
}
