show_MIAPE <- function(object) {
  if (any(c(length(object@instrumentModel) > 0,
            length(object@instrumentManufacturer) > 0,
            length(object@instrumentCustomisations) > 0))) {
    cat("Instrument : \n")
    cat("  Model:", countAndPrint(object@instrumentModel), "\n")
    cat("  Manufacturer:", countAndPrint(object@instrumentManufacturer), "\n")
    cat("  Customisations:", countAndPrint(object@instrumentCustomisations), "\n")
    cat("  Use 'msInfo(object)' for more MIAPE-MS information.\n")
  }
  ## code from original MIAPE show method
  tmp <- c("samples","preprocessing")
  Index <-c(length(object@samples) > 0,
            length(object@preprocessing) > 0)
  cat("Experiment data\n")
  cat("  Experimenter name:", object@name,"\n")
  cat("  Laboratory:", object@lab,"\n")
  cat("  Contact information:", object@contact,"\n")
  cat("  Contact email:", object@email,"\n")
  cat("  Title:", object@title,"\n")
  cat("  URL:", object@url,"\n")
  pmids <- pubMedIds(object)
  cat("  PMIDs:", pmids,"\n")
  if(length(object@abstract) > 0 && all(object@abstract!=""))
    cat("\n  Abstract: A", length(strsplit(object@abstract," ")[[1]]),
        "word abstract is available. Use 'abstract' method.\n")
  else
    cat("  No abstract available.\n")
  if (any(Index))
    cat("  Information is available on:", paste(tmp[Index],collapse=", "),"\n")
  nO <- notes(object)
  if (length(nO) > 0) {
    cat("  notes:\n" )
    if( is.list(nO) ) {
      nms <- names(nO)
      pw <- options("width")[[1]] - 6
      for(i in 1:length(nO) ) {
          cat("   ", nms[i], ":\n", sep="")
          if (all(is.na(nO[[i]]))) {
              cat("\n")
          } else {
              ni <- gsub("  ", "", gsub("\n", "", nO[[i]]))
              cat(strbreak(ni, width=pw, exdent=0), sep="\n      ")
          }
      }
    }
  }
}

msInfo_MIAPE <- function(object) {
  cat("MIAPE-MS information:\n")
  cat(" 1. General features:\n")
  cat("  Date stamp:", object@dateStamp,"\n")
  cat("  Contact:", object@contact,"\n")
  cat("  Name:", object@name,"\n")
  cat("  Laboratory:", object@lab,"\n")
  cat("  Instument model:", instrumentModel(object),"\n")
  cat("  Manufacturer:", instrumentManufacturer(object),"\n")
  cat("  Customisations:", instrumentCustomisations(object),"\n")
  cat("  Software:", object@softwareName,"\n")
  cat("  Version:", object@softwareVersion,"\n")
  cat("  Switching:", object@switchingCriteria,"\n")
  cat("  Param file:", object@parameterFile,"\n")
  cat(" 2. Ion source\n")
  cat("  Source:", ionSource(object),"\n")
  cat("  Source details:", ionSourceDetails(object),"\n")
  cat(" 3. Post-source componentry\n")
  cat("  Analyzer:", analyser(object),"\n")
  cat("  Analyzer details:", analyzerDetails(object),"\n")
  cat("  Collision gas:", object@collisionGas,"\n")
  cat("  Pressure:", object@collisionPressure," bars\n")
  cat("  Energy:", object@collisionEnergy,"\n")
  cat("  Detector type:", detectorType(object),"\n")
  cat("  Sensitivity:", object@detectorSensitivity,"\n")
  ## cat("See spectra for their description and\n")
  ## cat("the MSnProcess slot for processing details.\n")
}
