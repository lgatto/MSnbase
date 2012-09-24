show_MSnProcess <- function(object) {
  cat("- - - Processing information - - -\n")
  for (proc in object@processing)
    cat(proc,"\n")
  cat(" MSnbase version:",object@MSnbaseVersion,"\n")
}
