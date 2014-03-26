## mzTab - Reporting Proteomics Results
## ref: http://code.google.com/p/mztab/

mzTabDate <- function() format(Sys.time(), "%Y-%m-%dT%H%M%Z")

## Currently only supports export of unique units
  
writeCOM <- function(com, file = "")
  quite <- sapply(paste("COM", com), function(x) cat(x, file, "\n"))

##' \code{mzTab} is a light-weight, tab-delimited file format for proteomics data.
##' It describes general metadata, protein, peptide and small molecule information
##' (all of which are optinal), including quantitation and identification.
##' The metadata section (MTD) can be generated from an \code{\linkS4class{MSnSet}}
##' instance using \code{makeMTD}. The detailed description of all the parameters
##' can be found in the \code{mzTab} specification document (see references).
##' 
##' @title Creates the mzTab metadata section
##' @param x An instance of class \code{\linkS4class{MSnSet}}. 
##' @param unitId A \code{character} of lenth 1 or NULL (default),  in
##' which case \code{x}'s variable name will be used. This identifier
##' references the item under study all sections.
##' @param title A \code{character} of lenght 1 or \code{NULL}
##' (default), in which case \code{exptitle(x)} is used if available.
##' @param mtdDescription A \code{character} of length 1 describing the unit
##' or \code{NULL} (default) to ignore.
##' @param sampleProcessing A \code{list} of (possibly multiple) valid
##' \code{CVParam} objects or \code{NULL} (default) to ignore.
##' @param instrumentSource A list of valid \code{CVParam} instances
##' or \code{NULL} (default), in which case \code{ionSource(x)} is
##' used to generate a \code{CVParam}.
##' @param instrumentAnalyzer A list of valid \code{CVParam} instances
##' or \code{NULL} (default), in which case \code{analyzer(x)} is used
##' to generate a \code{CVParam}.
##' @param instrumentDetector A list of valid \code{CVParam} instances
##' or \code{NULL} (default), in which case \code{detectorType(x)} is
##' used to generate a \code{CVParam}.
##' @param software A list of valid \code{CVParam} instances describing
##' the ordered list of software used to process the data. 
##' \code{NULL} (default) to ignore.
##' @param fdr A list of valid \code{CVParam} instances describing the
##' unit's false discovery rate or \code{NULL} (default) to ignore.
##' @param publication A \code{character} (of lenght > 0) or
##' \code{NULL} (default), in which case \code{pubMedIds(x)} is used.
##' @param contactName A \code{character} (of length > 0) or
##' \code{NULL} (default), in which case \code{expinfo(x)["name"]} is
##' used.
##' @param contactAffiliation A \code{character} (of length > 0) or
##' \code{NULL} (default), in which case \code{expinfo(x)["lab"]} is
##' used. 
##' @param contactEmail A \code{character} (of length > 0) or
##' \code{NULL} (default), in which case \code{expemail(x)} is
##' used.
##' @param mtdUri A \code{character} (of length > 0) describing the unit's
##' uniform resource identifier (a PRIDE experiment or a PeptideAtlas build
##' for example). \code{NULL} (default) to ignore.
##' @param mtdModifications A list of (possibly multiplt) \code{CVParam} instances
##' describing all (distinct) PTMs reported in the unit. 
##' \code{NULL} (default) to ignore.
##' @param modProbabilityMethod A user defined \code{CVParam} reporting
##' the modifiction (position) probabilities. \code{NULL} (default) to ignore.
##' @param quantitationMethod A valid \code{CVParam}, a
##' \code{\linkS4class{ReporterIons}} instance or \code{NULL}
##' (default), in which case the isobaric tagging system is guessed
##' from the number of columns in \code{exprs(x)} (4 or 8 for iTRAQ, 6
##' for TMT).
##' @param protQuantUnit A valid \code{CVParam} or \code{NULL}
##' (default) to use PRIDE:0000330 (Arbitrary quantificatio unit).
##' @param pepQuantUnit A valid \code{CVParam} or \code{NULL}
##' (default) to use PRIDE:0000330 (Arbitrary quantificatio unit).
##' @param msFileFormat A list of valid \code{CVParam} instances to
##' \code{NULL} (default), in which case, the extension of
##' \code{fileNames(x)[1]} is used to define the appropriate
##' \code{CVParam}. Recognised extensions are \code{mzData},
##' \code{mzXML}, \code{mzML} or \code{mgf}.
##' @param msFileLocation A \code{character} (of length > 0) or
##' \code{NULL} (default), in which case \code{fileNames(x)} is used.
##' @param msFileIdFormat A list of \code{CVParam} instances describing
##' the original identification format used in the external data file.
##' \code{NULL} (default) to ignore.
##' @param custom A list of user defined \code{CVParam} instances with
##' additional parameters describing the unit. \code{NULL} (default) to ignore.
##' @param species_ A list of (possibly several) \code{CVParam}
##' instances with the respective (sub-)unit species.
##' \code{NULL} (default) to ignore.
##' @param tissue_ A list of (possibly several) \code{CVParam}
##' instances describing the respective (sub-)unit tissue.
##' \code{NULL} (default) to ignore.
##' @param cellType_ A list of (possibly several) \code{CVParam}
##' instances describing the respective (sub-)unit cell type.
##' \code{NULL} (default) to ignore.
##' @param disease_ A list of (possibly several) \code{CVParam}
##' instances describing the respective (sub-)unit disease states.
##' \code{NULL} (default) to ignore.
##' @param description_ A list of \code{characters} describing the (sub-)unit
##' in human raedable free text. \code{NULL} (default) to ignore.
##' @param quantitationReagent_ A list of \code{CVParam} instances or
##' \code{NULL} (default), in which case the reporter ions as defined
##' by \code{quantitationMethod} as used. 
##' @param custom_ A list of user defined \code{CVParam} instances with
##' additional (sub-)unit properties. \code{NULL} (default) to ignore. 
##' @return A \code{character} defining the mzTab metadata section.
##' @references mzTab - Reporting Proteomics Results (\url{http://code.google.com/p/mztab/})
##' @seealso \code{\link{makePEP}} and \code{\link{makePRT}} to generate
##' \code{mzTab} peptide and protein sections.
##' @author Laurent Gatto
makeMTD <- function(x,
                    unitId = NULL,
                    title = NULL,
                    mtdDescription = NULL,
                    sampleProcessing = NULL,
                    instrumentSource = NULL,
                    instrumentAnalyzer = NULL,
                    instrumentDetector = NULL,
                    software = NULL,
                    fdr = NULL,     
                    publication = NULL,
                    contactName = NULL,
                    contactAffiliation = NULL,
                    contactEmail = NULL,
                    mtdUri = NULL,
                    mtdModifications = NULL,
                    modProbabilityMethod = NULL,
                    quantitationMethod = NULL,  
                    protQuantUnit = NULL, 
                    pepQuantUnit = NULL,  
                    msFileFormat = NULL,  
                    msFileLocation = NULL,
                    msFileIdFormat = NULL,
                    custom = NULL,
                    species_ = NULL,
                    tissue_ = NULL,
                    cellType_ = NULL,
                    disease_ = NULL,  
                    description_ = NULL,
                    quantitationReagent_ = NULL, 
                    custom_  = NULL 
                    ) {
  if (class(x) != "MSnSet")
    stop("Need instance of class 'MSnSet' to write to mzTab format.")
  
  addToMtd <- function(x, k,
                       .unitId = unitId,
                       .mtd = mtd,
                       pre = "MTD", sep = "\t") {
    if (is.null(x)) {
      ans <- .mtd
    } else {      
      k <- paste0(.unitId, "-", k)
      x <- as.character(x) ## to convert CVParam's 
      x <- paste(pre, k, x, sep = sep)
      x <- paste0(x, "\n")
      ans <- append(.mtd, x)
    }    
    return(ans)
  }

  .n <- 0 ## number of reporter ions
  mtd <- c()
  if (is.null(unitId)) 
    unitId <- getVariableName(match.call(), "x")
  if (length(unitId) > 1) {
    unitId <- unitId[1]
    warning("Expecting 'unitId' of length 1. Dropping others.")
  }
  
  if (is.null(title)) {
    if (experimentData(x)@title != "") 
      title <- experimentData(x)@title
  }
  mtd <- addToMtd(title, "title")  
  mtd <- addToMtd(mtdDescription, "description")

  mtd <- addToMtd(sampleProcessing,
                  paste0("sample_processing[", 1:length(sampleProcessing), "]")) 

  if (is.null(instrumentSource)) {
    if (length(ionSource(x)) > 0) {
      instrumentSource <- ionSource(x)
      instrumentSource <- olsQuery(instrumentSource, "PSI", exact = TRUE)
      if (length(instrumentSource) != 1) {
        instrumentSource <- NULL
        warning("No unique instrumentSource ontology term found.")
      } else {
        .label <- strsplit(names(instrumentSource), ":")[[1]][1]
        .accession <- names(instrumentSource)
        .name <- term(.accession, .label)
        instrumentSource <- new("CVParam", label = .label, accession = .accession, name = .name)
      }
    }
  }
  mtd <- addToMtd(as.character(instrumentSource), "instrument[1]-source")

  if (is.null(instrumentAnalyzer)) {
    if (length(analyzer(x)) > 0) {
      instrumentAnalyzer <- analyzer(x)
      instrumentAnalyzer <- olsQuery(instrumentAnalyzer, "MS", exact = TRUE)
      if (length(instrumentAnalyzer) != 1) {
        instrumentAnalyzer <- NULL
        warning("No unique instrumentAnalyzer ontology term found.")
      } else {
        .label <- strsplit(names(instrumentAnalyzer), ":")[[1]][1]
        .accession <- names(instrumentAnalyzer)
        .name <- term(.accession, .label)
        instrumentAnalyzer <- new("CVParam",
                                  label = .label,
                                  accession = .accession,
                                  name = .name)
      }
    }
  }
  mtd <- addToMtd(as.character(instrumentAnalyzer), "instrument[1]-analyzer")

  if (is.null(instrumentDetector)) {
    if (length(detectorType(x)) > 0) {
      instrumentDetector <- detectorType(x)
      instrumentDetector <- olsQuery(instrumentDetector, "MS", exact = TRUE)
      if (length(instrumentDetector) != 1) {
        instrumentDetector <- NULL
        warning("No unique instrumentDetector ontology term found.")
      } else {
        .label <- strsplit(names(instrumentDetector), ":")[[1]][1]
        .accession <- names(instrumentDetector)
        .name <- term(.accession, .label)
        instrumentDetector <- new("CVParam", label = .label, accession = .accession, name = .name)
      }
    }
  }
  mtd <- addToMtd(as.character(instrumentDetector), "instrument[1]-detector")
  
  mtd <- addToMtd(software, paste0("software[", 1:length(software), "]")) 

  mtd <- addToMtd(fdr, paste0("fdr[", 1:length(fdr), "]"))
  
  if (is.null(publication)) {
    if (experimentData(x)@pubMedIds != "") {
      publication <- pubMedIds(x)
      publication <- paste0("pubmed:", publication)
    }
  }
  mtd <- addToMtd(publication, paste0("publication[", 1:length(publication), "]"))

  if (is.null(contactName)) {
    if (expinfo(x)["name"] != "") 
      contactName <- expinfo(x)["name"]
  }
  mtd <- addToMtd(contactName, paste0("contact[", 1:length(contactName), "]-name"))

  if (is.null(contactAffiliation)) {
    if (expinfo(x)["lab"] != "") 
      contactAffiliation <- expinfo(x)["lab"]
  }
  mtd <- addToMtd(contactAffiliation,
                  paste0("contact[", 1:length(contactAffiliation), "]-affiliation"))

  if (is.null(contactEmail)) {
    if (expemail(x) != "") 
      contactEmail <- expemail(x)
  }
  mtd <- addToMtd(contactEmail,
                  paste0("contact[", 1:length(contactEmail), "]-email"))

  mtd <- addToMtd(mtdUri, paste0("uri[", 1:length(mtdUri), "]"))
  mtd <- addToMtd(paste(sapply(mtdModifications, as.character), collapse = " | "),
                  "mod")
  mtd <- addToMtd(modProbabilityMethod, "mod-probability_method")

  if (is.null(quantitationMethod)) {
    ## Guessing the reporters from ncol
    if (ncol(x) == 6) {
      quantitationMethod <- new("CVParam", label = "PRIDE", accession = "PRIDE:0000314", name = "TMT")
      .n <- 6
    }
    if (ncol(x) %in% c(4, 8)) {
      quantitationMethod <- new("CVParam", label = "PRIDE", accession = "PRIDE:0000313", name = "iTRAQ")
      .n <- ncol(x)
    }      
  } else if (class(quantitationMethod) == "ReporterIons") {
    .n <- length(quantitationMethod)
    quantitationMethod <- substr(names(quantitationMethod),
                                 1, nchar(names(quantitationMethod)) - 1)    
    if (quantitationMethod == "iTRAQ") {
      quantitationMethod <- new("CVParam", label = "PRIDE", accession = "PRIDE:0000313", name = "iTRAQ")
    } else {
      quantitationMethod <- new("CVParam", label = "PRIDE", accession = "PRIDE:0000314", name = "TMT")
    }
  }
  if (!is.null(quantitationMethod)) ## then is must be a CVParam
    mtd <- addToMtd(as.character(quantitationMethod), "quantitation_method")

  ## default quant unit is PRIDE:0000330, Arbitrary quantificatio unit
  if (is.null(protQuantUnit)) {
    .accession <- "PRIDE:0000330"
    .name <- term(.accession, "PRIDE")
    protQuantUnit <- new("CVParam", label = "PRIDE", accession = .accession, name = .name)
  }
  mtd <- addToMtd(as.character(protQuantUnit), "protein-quantification_unit")

  if (is.null(pepQuantUnit)) {
    .accession <- "PRIDE:0000330"
    .name <- term(.accession, "PRIDE")
    pepQuantUnit <- new("CVParam", label = "PRIDE", accession = .accession, name = .name)    
  }
  mtd <- addToMtd(as.character(pepQuantUnit), "peptide-quantification_unit")

  if (is.null(msFileFormat)) {
    if (length(fileNames(x)) > 0) {
      ext <- tail(unlist(strsplit(fileNames(x), "\\.")), n = 1)
      if (tolower(ext) %in% c("mzml", "mzxml", "mzdata", "mgf")) {
        .name <- olsQuery(paste(ext, "format"), "MS")
        .accession <- names(.name)
        msFileFormat <- new("CVParam", label = "MS", accession = .accession, name = .name)
      ## } else if (tolower(ext) == "mgf") {
      ##   .name <- olsQuery("peak list scans", "MS", exact = TRUE)
      ##   .accession <- names(.name)
      ##   msFileFormat <- new("CVParam", label = "MS", accession = .accession, name = .nam)e
      } else {
        warning("File format '", ext, "' not recognised.")
        msFileFormat <- NULL
      }
    }
  }
  if (!is.null(msFileFormat))
    mtd <- addToMtd(as.character(msFileFormat), paste0("ms_file[", 1:length(ext), "]-format"))
  
  if (is.null(msFileLocation)) {
    if (length(fileNames(x)) > 0) 
      msFileLocation <- fileNames(x)
  }
  mtd <- addToMtd(msFileLocation,
                  paste0("ms_file[", 1:length(msFileLocation), "]-location"))
  mtd <- addToMtd(msFileIdFormat,
                  paste0("ms_file[", 1:length(msFileIdFormat), "]-id_format"))
  mtd <- addToMtd(custom,
                  paste0("custom[", 1:length(custom), "]"))

  if (!is.null(species_)) {
    for (i in 1:.n) {
      sub <- species_[[i]]
      if (length(sub) == 1) {
        mtd <- addToMtd(as.character(sub),
                        paste0("sub[", i, "]-species[1]"))
      } else {
        mtd <- addToMtd(sapply(sub, as.character),
                        paste0("sub[", i, "]-species[", 1:length(sub), "]"))
      }
    }
  }

  if (!is.null(tissue_)) {
    for (i in 1:.n) {
      sub <- tissue_[[i]]
      if (length(sub) == 1) {
        mtd <- addToMtd(as.character(sub),
                        paste0("sub[", i, "]-tissue[1]"))
      } else {
        mtd <- addToMtd(sapply(sub, as.character),
                        paste0("sub[", i, "]-tissue[", 1:length(sub), "]"))
      }
    }
  }

  if (!is.null(cellType_)) {
    for (i in 1:.n) {
      sub <- cellType_[[i]]
      if (length(sub) == 1) {
        mtd <- addToMtd(as.character(sub),
                        paste0("sub[", i, "]-cell_type[1]"))
      } else {
        mtd <- addToMtd(sapply(sub, as.character),
                        paste0("sub[", i, "]-cell_type[", 1:length(sub), "]"))
      }
    }
  }
  
  if (!is.null(disease_)) {
    for (i in 1:.n) {
      sub <- disease_[[i]]
      if (length(sub) == 1) {
        mtd <- addToMtd(as.character(sub),
                        paste0("sub[", i, "]-disease[1]"))
      } else {
        mtd <- addToMtd(sapply(sub, as.character),
                        paste0("sub[", i, "]-disease[", 1:length(sub), "]"))
      }
    }
  }

  if (!is.null(description_))
    mtd <- addToMtd(sapply(description_, as.character),
                    paste0("sub[", 1:length(description_), "]-description"))

  if (!is.null(custom_))
    mtd <- addToMtd(sapply(custom_, as.character),
                    paste0("sub[", 1:length(custom_), "]-custom")) 

  if (!is.null(quantitationMethod) & is.null(quantitationReagent_)) {
    .accession <- .names <- c()    
    if (.n == 4) {
      .accessions <- c("PRIDE:0000264", "PRIDE:0000114", "PRIDE:0000115", "PRIDE:0000116")
      .names <- c("iTRAQ reagent 113", "iTRAQ reagent 114", "iTRAQ reagent 115", "iTRAQ reagent 116")
      quantitationReagent_ <- TRUE
    } else if (.n == 8) {
      .accessions <- c("PRIDE:0000264", "PRIDE:0000114", "PRIDE:0000115", "PRIDE:0000116",
                       "PRIDE:0000117", "PRIDE:0000265", "PRIDE:0000266", "PRIDE:0000167")
      .names <- c("iTRAQ reagent 113", "iTRAQ reagent 114", "iTRAQ reagent 115", "iTRAQ reagent 116",
                  "iTRAQ reagent 117", "iTRAQ reagent 118", "iTRAQ reagent 119", "iTRAQ reagent 121")
      quantitationReagent_ <- TRUE
    } else {  ## .n == 8) 
      .accessions <- c("PRIDE:0000285", "PRIDE:0000286", "PRIDE:0000287",
                       "PRIDE:0000288", "PRIDE:0000289", "PRIDE:0000290")
      .names <- c("TMT reagent 126", "TMT reagent 127", "TMT reagent 128",
                  "TMT reagent 129", "TMT reagent 130", "TMT reagent 131")
      quantitationReagent_ <- TRUE
    }
    if (isTRUE(quantitationReagent_)) {
      quantitationReagent_<- mapply(function(a, b)
                                    new("CVParam", label = "PRIDE", accession = a, name = b),
                                    .accessions, .names)
    }

    if (is.list(quantitationReagent_)) {
      for (i in 1:.n)
        mtd <- addToMtd(as.character(quantitationReagent_[[i]]),
                        paste0("sub[", i, "]-quantitation_reagent"))
    }

  }
  return(mtd)
}



##' \code{mzTab} is a light-weight, tab-delimited file format for proteomics data.
##' It describes general metadata, protein, peptide and small molecule information
##' (all of which are optinal), including quantitation and identification.
##' The peptide section (PEH header and PEP tabular data) can be generated from
##' an \code{\linkS4class{MSnSet}}  instance using \code{makePEP}.
##' The detailed description of all the parameters can be found in the \code{mzTab}
##' specification document (see references).
##'
##' @title Creates the mzTab peptide section
##' @param x An instance of class \code{\linkS4class{MSnSet}}. 
##' @param sequence A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with the
##' peptide sequence. Default is \code{NA}.
##' @param pepAccession A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with the assigned
##' protein accession. Default is \code{NA}.
##' @param unitId A \code{character} of lenth 1 or NULL (default), in
##' which case \code{x}'s variable name will be used.
##' @param unique A \code{logical} (converted to numeric to comply
##' with format specification) of \code{length(nrow(x)} (will be
##' recycled a whole number of times if of different length) specifying
##' if peptide is proteotypic. Default is \code{NA}.
##' @param pepDatabase A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) describing
##' the protein database used for peptide identification.
##' Default is \code{NA}.
##' @param pepDatabaseVersion A \code{character} of length \code{nrow(x)}
##' (will be recycled a whole number of times if of different length) with
##' the database version. Default is \code{NA}.
##' @param pepSearchEngine A \code{list} of length \code{nrow(x)} (of possibly
##' multiple lists of) \code{CVParam} instances identifying the search engine
##' used for peptide identification. Default is \code{NA}.
##' @param pepSearchEngineScore A \code{list} of length \code{nrow(x)} (of possibly
##' multiple lists of) \code{CVParam} instances specifying peptide identification
##' scores. Default is \code{NA}.
##' @param pepReliability A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length). Values should
##' be \code{1} (high reliability), \code{2} (medium reliability) or \code{3}
##' (poor reliability). Default is \code{NA}.
##' @param pepModifications A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) describing the
##' modifications and their position (see \code{mzTab} format specifications
##' for details). Default is \code{NA}.
##' @param retentionTime A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length). Note that currently,
##' unique retention times are expected, but could be extended to multiple
##' times. Default is \code{NA}.
##' @param charge A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) indicating
##' peptide charge state.  Default is \code{NA}.
##' @param massToCharge A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with the peptides
##' precursor mass to charge ratio. Default is \code{NA}.
##' @param pepUri A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with
##' peptide uniform resource identifiers (link to PRIDE database for instance). 
##' Default is \code{NA}.
##' @param spectraRef A \code{character} in the format \code{ms_file[1-n]:{SPEC_REF}}
##' (see \code{mzTab} specifications for details) of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length).
##' Default is \code{NA}.
##' @param pepAbundance A \code{numeric} of length \code{nrow(x}) or
##' \code{matrix} with \code{nrow(x)} rows if multiple sub-samples are
##' reported (see metadata section), specifying the peptides abundance.
##' If \code{NULL} (default), ignored.
##' @param pepAbundanceStdev A \code{numeric} of length \code{nrow(x}) or
##' \code{matrix} with \code{nrow(x)} rows if multiple sub-samples are
##' reported (see metadata section), specifying the standard deviation of
##' peptides abundances. If \code{NULL} (default), ignored.
##' If \code{pepAbundance} is not \code{NULL}, then \code{pepAbundanceStdev}
##' is \code{NA} if not specified.
##' @param pepAbundanceSterr A \code{numeric} of length \code{nrow(x}) or
##' \code{matrix} with \code{nrow(x)} rows if multiple sub-samples are
##' reported (see metadata section), specifying the standard error of
##' peptides abundances. If \code{NULL} (default), ignored.
##' If \code{pepAbundance} is not \code{NULL}, then \code{pepAbundanceSterr}
##' is \code{NA} if not specified.
##' @param pepOpt_ An optional \code{character} of \code{character}
##' \code{matrix} (possibly populated with text represenations of
##' \code{CVParam} instances) for any custom peptide annotation.
##' Default is \code{NULL} to ignore.
##' @return A \code{data.frame} defining the \code{mzTab} peptide section.
##' @seealso \code{\link{makeMTD}} and \code{\link{makePRT}} to generate
##' \code{mzTab} metadata and protein sections.
##' @author Laurent Gatto
makePEP <- function(x,
                    sequence = NA,
                    pepAccession = NA,
                    unitId = NULL,
                    unique = NA,
                    pepDatabase = NA,
                    pepDatabaseVersion = NA,
                    pepSearchEngine = NA,
                    pepSearchEngineScore = NA,
                    pepReliability = NA,
                    pepModifications = NA,
                    retentionTime = NA,
                    charge = NA,
                    massToCharge = NA,
                    pepUri = NA,
                    spectraRef = NA,
                    pepAbundance = NULL,
                    pepAbundanceStdev = NULL,
                    pepAbundanceSterr = NULL,
                    pepOpt_ = NULL) {
  if (class(x) != "MSnSet")
    stop("Need instance of class 'MSnSet' to write to mzTab format.")  
  
  if (is.null(unitId)) 
    unitId <- getVariableName(match.call(), "x")

  if (!all(is.na(pepSearchEngine))) {
    pepSearchEngine <- sapply(pepSearchEngine,
                           function(.seng) {
                             if (is.list(.seng)) {
                               .classes <- sapply(.seng, class)
                               if (any(sapply(.seng, class) != "CVParam")) {
                                 .classes <- .classes[.classes != "CVParam"]
                                 .classes <- paste(.classes, collapse = ", ")
                                 stop("Expecting list of 'CVParam' instances in 'pepSearchEngine'.",
                                      " Found ", .classes, ".")
                               }
                               return(paste(sapply(.seng, as.character), collapse = " |"))
                             } else {
                               if (class(.seng) != "CVParam")
                                 stop("Expecting list of 'CVParam' instances in 'pepSearchEngine'.",
                                      " Found ", class(.seng), ".")
                               return(as.character(.seng))
                             }
                           })
  }

  if (!all(is.na(pepSearchEngineScore))) {
    pepSearchEngineScore <- sapply(pepSearchEngineScore,
                           function(.sscore) {
                             if (is.list(.sscore)) {
                               .classes <- sapply(.sscore, class)
                               if (any(sapply(.sscore, class) != "CVParam")) {
                                 .classes <- .classes[.classes != "CVParam"]
                                 .classes <- paste(.classes, collapse = ", ")
                                 stop("Expecting list of 'CVParam' instances in 'pepSearchEngineScore'.",
                                      " Found ", .classes, ".")
                               }
                               return(paste(sapply(.sscore, as.character), collapse = " |"))
                             } else {
                               if (class(.sscore) != "CVParam")
                                 stop("Expecting list of 'CVParam' instances in 'pepSearchEngineScore'.",
                                      " Found ", class(.sscore), ".")
                               return(as.character(.sscore))
                             }
                           })
  }

  .__pep <- data.frame(sequence = sequence,
                      accession = pepAccession,
                      unit_id = unitId,
                      unique = unique,
                      database = pepDatabase,
                      database_version = pepDatabaseVersion,
                      search_engine = pepSearchEngine,
                      search_engine_score = pepSearchEngineScore,
                      reliability = pepReliability,
                      modifications = pepModifications,
                      retention_time = retentionTime,
                      charge = charge,
                      mass_to_charge = massToCharge,
                      uri = pepUri,
                      spectra_ref = spectraRef)

  if (!is.null(pepAbundance)) {
    if (!is.matrix(pepAbundance)) 
      pepAbundance <- matrix(pepAbundance, ncol = 1)
    colnames(pepAbundance) <- paste0("peptide_abundance_sub[", 1:ncol(pepAbundance), "]")
    
    if (is.null(pepAbundanceStdev)) {
      pepAbundanceStdev <- matrix(NA, nrow = nrow(pepAbundance), ncol = ncol(pepAbundance))
    } else {
      if (!is.matrix(pepAbundanceStdev))
        pepAbundanceStdev <- matrix(pepAbundanceStdev, ncol = 1)
    }
    colnames(pepAbundanceStdev) <- paste0("peptide_abundance_stdev_sub[", 1:ncol(pepAbundanceStdev), "]") 

    if (is.null(pepAbundanceSterr)) {
      pepAbundanceSterr <- matrix(NA, nrow = nrow(pepAbundance), ncol = ncol(pepAbundance))
    } else {
      if (!is.matrix(pepAbundanceSterr))
        pepAbundanceSterr <- matrix(pepAbundanceSterr, ncol = 1)
    }
    colnames(pepAbundanceSterr) <- paste0("peptide_abundance_std_error_sub[", 1:ncol(pepAbundanceSterr), "]")      

    
    if (ncol(pepAbundanceStdev) != ncol(pepAbundance))
      stop("pepAbundance and pepAbundanceStdev have different number of columns.")
    if (ncol(pepAbundanceSterr) != ncol(pepAbundance))
      stop("pepAbundance and pepAbundanceSterr have different number of columns.")

    .__pep <- cbind(.__pep, pepAbundance, pepAbundanceStdev, pepAbundanceSterr)
  }

  if (!is.null(pepOpt_)) {
    if (!all(subset(colnames(pepOpt_), 1, 4) == "opt_"))
      colnames(pepOpt_) <- paste0("opt_", colnames(pepOpt_))
    .__pep <- cbind(.__pep, pepOpt_)
  }
  
  return(.__pep)
}


##' \code{mzTab} is a light-weight, tab-delimited file format for proteomics data.
##' It describes general metadata, protein, peptide and small molecule information
##' (all of which are optinal), including quantitation and identification.
##' The proteine section (PRH header and PRT tabular data) can be generated from
##' an \code{\linkS4class{MSnSet}}  instance using \code{makePRT}.
##' The detailed description of all the parameters can be found in the \code{mzTab}
##' specification document (see references).
##'
##' @title Creates the mzTab protein section
##' @param x An instance of class \code{\linkS4class{MSnSet}}. 
##' @param protAccession A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with the 
##' protein accession. Default is \code{NA}.
##' @param unitId A \code{character} of lenth 1 or NULL (default),  in
##' which case \code{x}'s variable name will be used.
##' @param protDescription A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with the 
##' protein name or description. Default is \code{NA}.
##' @param taxId A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) referencing
##' the species NCBI/NEWT taxonomy id. Default is \code{NA}.
##' @param species A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) describing the
##' species in human readable form. Default is \code{NA}.
##' @param protDatabase A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) describing the
##' protein database. Default is \code{NA}.
##' @param protDatabaseVersion A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) describing the
##' database version. Default is \code{NA}.
##' @param protSearchEngine A \code{list} of length \code{nrow(x)}
##' (of possibly several) \code{CVParam} instances describing the search engine
##' used for protein identification. Default is \code{NA}.
##' @param protSearchEngineScore A \code{list} of length \code{nrow(x)} (of possibly
##' multiple lists of) \code{CVParam} instances specifying peptide identification
##' scores. Default is \code{NA}.
##' @param protReliability A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length). Values should
##' be \code{1} (high reliability), \code{2} (medium reliability) or \code{3}
##' (poor reliability). Default is \code{NA}.
##' @param numPep A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) indicating the number
##' of peptides identifying the proteins. Default is \code{NA}.
##' @param numPepDistinct A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) indicating the number
##' of distinct peptides (sequence and modifications) identifying the proteins.
##' Default is \code{NA}.
##' @param numPepUnambiguous A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) indicating the number
##' of unambiguous disctinct peptides identifying the proteins. Default is \code{NA}.
##' @param ambiguityMembers A \code{character} of comma-separated protein accessions. 
##' See the \code{mzTab} specification document for details. Defaut is \code{NA}.
##' @param protModifications A \code{character} of comma-delimited modifications/scores/positions 
##' describing the proteins. See the \code{mzTab}specification document for details.
##' Defaut is \code{NA}.
##' @param protUri A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with
##' peptide uniform resource identifiers (link to PRIDE database for instance). 
##' Default is \code{NA}.
##' @param goTerms A \code{character} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with comma-delimited
##' GO terms describing the proteins. Default is \code{NA}.
##' @param protCoverage A \code{numeric} of length \code{nrow(x)} (will be
##' recycled a whole number of times if of different length) with the protein
##' coverages ranging between 0 and 1. Default is \code{NA}.
##' @param protAbundance A \code{numeric} of length \code{nrow(x}) or
##' \code{matrix} with \code{nrow(x)} rows if multiple sub-samples are
##' reported (see metadata section), specifying the protein abundance.
##' If \code{NULL} (default), ignored.
##' @param protAbundanceStdev A \code{numeric} of length \code{nrow(x}) or
##' \code{matrix} with \code{nrow(x)} rows if multiple sub-samples are
##' reported (see metadata section), specifying the standard deviation of
##' protein abundances. If \code{NULL} (default), ignored.
##' If \code{protAbundance} is not \code{NULL}, then \code{protAbundanceStdev}
##' is \code{NA} if not specified.
##' @param protAbundanceSterr A \code{numeric} of length \code{nrow(x}) or
##' \code{matrix} with \code{nrow(x)} rows if multiple sub-samples are
##' reported (see metadata section), specifying the standard error of
##' protein abundances. If \code{NULL} (default), ignored.
##' If \code{protAbundance} is not \code{NULL}, then \code{protAbundanceSterr}
##' is \code{NA} if not specified.
##' @param protOpt_ An optional \code{character} of \code{character}
##' \code{matrix} (possibly populated with text represenations of
##' \code{CVParam} instances) for any custom protein annotation.
##' Default is \code{NULL} to ignore.
##' @return A \code{data.frame} defining the \code{mzTab} protein section.
##' @references mzTab - Reporting Proteomics Results (\url{http://code.google.com/p/mztab/})
##' @seealso \code{\link{makeMTD}} and \code{\link{makePEP}} to generate
##' \code{mzTab} metadat and peptide sections.
##' @author Laurent Gatto
makePRT <- function(x,
                    protAccession = NA,
                    unitId = NULL,
                    protDescription = NA,
                    taxId = NA,
                    species = NA,
                    protDatabase = NA,
                    protDatabaseVersion = NA,
                    protSearchEngine = NA,
                    protSearchEngineScore = NA,
                    protReliability = NA,
                    numPep = NA,
                    numPepDistinct = NA,
                    numPepUnambiguous = NA,
                    ambiguityMembers = NA,
                    protModifications = NA,
                    protUri = NA,
                    goTerms = NA,
                    protCoverage = NA,
                    protAbundance = NULL,
                    protAbundanceStdev = NULL,
                    protAbundanceSterr = NULL,
                    protOpt_ = NULL) {
  if (is.null(unitId)) 
    unitId <- getVariableName(match.call(), "x")
  
  if (!all(is.na(protSearchEngine))) {
    protSearchEngine <- sapply(protSearchEngine,
                               function(.seng) {
                                 if (is.list(.seng)) {
                                   .classes <- sapply(.seng, class)
                                   if (any(sapply(.seng, class) != "CVParam")) {
                                     .classes <- .classes[.classes != "CVParam"]
                                     .classes <- paste(.classes, collapse = ", ")
                                     stop("Expecting list of 'CVParam' instances in 'protSearchEngine'.",
                                          " Found ", .classes, ".")
                                   }
                               return(paste(sapply(.seng, as.character), collapse = " |"))
                             } else {
                               if (class(.seng) != "CVParam")
                                 stop("Expecting list of 'CVParam' instances in 'protSearchEngine'.",
                                      " Found ", class(.seng), ".")
                               return(as.character(.seng))
                             }
                           })
  }

  if (!all(is.na(protSearchEngineScore))) {
    protSearchEngineScore <- sapply(protSearchEngineScore,
                           function(.sscore) {
                             if (is.list(.sscore)) {
                               .classes <- sapply(.sscore, class)
                               if (any(sapply(.sscore, class) != "CVParam")) {
                                 .classes <- .classes[.classes != "CVParam"]
                                 .classes <- paste(.classes, collapse = ", ")
                                 stop("Expecting list of 'CVParam' instances in 'protSearchEngineScore'.",
                                      " Found ", .classes, ".")
                               }
                               return(paste(sapply(.sscore, as.character), collapse = " |"))
                             } else {
                               if (class(.sscore) != "CVParam")
                                 stop("Expecting list of 'CVParam' instances in 'protSearchEngineScore'.",
                                      " Found ", class(.sscore), ".")
                               return(as.character(.sscore))
                             }
                           })
  }
  
  .__prot <- data.frame(accession = protAccession,
                        unit_id = unitId,
                        description = protDescription,
                        taxid = taxId,
                        species = species,
                        database = protDatabase,
                        database_version = protDatabaseVersion,
                        search_engine = protSearchEngine,
                        search_engine_score = protSearchEngineScore,
                        reliability = protReliability,
                        num_peptides = numPep,
                        num_peptides_distinct = numPepDistinct,
                        num_peptides_unambiguous = numPepUnambiguous,
                        ambiguity_members = ambiguityMembers,
                        modifications = protModifications,
                        uri = protUri,
                        go_terms = goTerms,
                        protein_coverage = protCoverage)
                     
  
  if (!is.null(protAbundance)) {
    if (!is.matrix(protAbundance)) 
      protAbundance <- matrix(protAbundance, ncol = 1)
    colnames(protAbundance) <- paste0("protein_abundance_sub[", 1:ncol(protAbundance), "]")
    
    if (is.null(protAbundanceStdev)) {
      protAbundanceStdev <- matrix(NA, nrow = nrow(protAbundance), ncol = ncol(protAbundance))
    } else {
      if (!is.matrix(protAbundanceStdev))
        protAbundanceStdev <- matrix(protAbundanceStdev, ncol = 1)
    }
    colnames(protAbundanceStdev) <- paste0("protein_abundance_stdev_sub[", 1:ncol(protAbundanceStdev), "]") 

    if (is.null(protAbundanceSterr)) {
      protAbundanceSterr <- matrix(NA, nrow = nrow(protAbundance), ncol = ncol(protAbundance))
    } else {
      if (!is.matrix(protAbundanceSterr))
        protAbundanceSterr <- matrix(protAbundanceSterr, ncol = 1)
    }
    colnames(protAbundanceSterr) <- paste0("protein_abundance_std_error_sub[", 1:ncol(protAbundanceSterr), "]")      

    
    if (ncol(protAbundanceStdev) != ncol(protAbundance))
      stop("protAbundance and protAbundanceStdev have different number of columns.")
    if (ncol(protAbundanceSterr) != ncol(protAbundance))
      stop("protAbundance and protAbundanceSterr have different number of columns.")

    .__prot <- cbind(.__prot, protAbundance, protAbundanceStdev, protAbundanceSterr)
  }

  if (!is.null(protOpt_)) {
    if (!all(subset(colnames(protOpt_), 1, 4) == "opt_"))
      colnames(protOpt_) <- paste0("opt_", colnames(protOpt_))
    .__prot <- cbind(.__prot, protOpt_)
  }
  
  return(.__prot)
}


##' This function generates an \code{mzTab} file based on the data
##' available in teh \code{x} \code{MSnSet} instance and additional 
##' information passed by the user. It make used of the respective
##' section generators to create appropriate metadata, peptide and
##' protein sections. If peptide and protein sections need to be
##' generated, one has to first create the \code{mzTab} file with
##' (metadata, optional, but recommended and) protein data and then
##' append the peptide data, as per \code{mzTab} specification
##' (\url{http://code.google.com/p/mztab/}). See the example section.
##'
##' @title Writes an 'MSnSet' to an mzTab file
##' @param x An instance of class \code{MSnSet}.
##' @param what One of \code{"PEP"} or \code{"PRT"} defining whether
##' peptide or protein data is to be saved.
##' @param append Logical. Should the data be appended to \code{file}.
##' Default is \code{FALSE}.
##' @param MTD Logical. Should the metadata section be generated.
##' Default is \code{TRUE}.
##' @param file A \code{character} naming the file to print to.
##' @param ... Additional parameters passed to the respective section
##' generators: \code{\link{makeMTD}} for metadata, \code{\link{makePEP}}
##' for peptides and \code{\link{makePRT}} for proteins.
##' @return None (invisible \code{NULL}).
##' @author Laurent Gatto
##' @seealso Functions to generate metadata (\code{\link{makeMTD}}),
##' peptide data (\code{\link{makePEP}}) and proteins
##' (\code{\link{makePRT}}). \code{\link{readMzTabData}} to create
##' \code{"\linkS4class{MSnSet}"} instances from an \code{mzTab} file.
##' @references The \code{mzTab} specification document and example
##' files: \url{http://code.google.com/p/mztab/}.
##' @examples
##' mzTabFile <- tempfile()
##' data(itraqdata)
##' pep <- quantify(itraqdata, reporters = iTRAQ4)
##' prot <- combineFeatures(pep, groupBy = fData(pep)$ProteinAccession)
##' fvarLabels(pep)
##' ## First write metadata and protein data
##' writeMzTabData(prot, what = "PRT", file = mzTabFile,
##'                append = FALSE, MTD = TRUE,
##'                protAccession = fData(prot)$ProteinAccession,
##'                protDescription = fData(prot)$ProteinDescription,
##'                protAbundance = exprs(prot))
##' ## append peptide data, without metadata section
##' writeMzTabData(pep, what = "PEP", file = mzTabFile,
##'                append = TRUE, MTD = FALSE,
##'                sequence = fData(pep)$PeptideSequence,
##'                charge = fData(pep)$charge,
##'                retentionTime = fData(pep)$retention.time,
##'                pepAbundance = exprs(pep))
writeMzTabData <- function(x,
                           what = c("PEP", "PRT"),
                           append = FALSE,
                           MTD = TRUE,
                           file, ...) {
    warning("Support for mzTab version 0.9 only. Support will be added soon.")
    if (!require(rols))
        stop("The 'rols' package is required for mzTab write support.")
    
    if (missing(file))
        stop("To which file would you like the data to be saved to?")
    
    what <- match.arg(what)

    params <- as.list(match.call())[-1]
    params <- params[!names(params) %in% c("what", "file", "MTD", "append")]
    if (is.null(params$unitId))
        params$unitId <- getVariableName(match.call(), "x")    
    paramNames <- names(params)

    selMTD <- paramNames %in% names(formals(makeMTD))
    selPEP <- paramNames %in% names(formals(makePEP))
    selPRT <- paramNames %in% names(formals(makePRT))
    
    noSel <- !(selMTD | selPEP | selPRT)
    if (any(noSel))
        warning("Parameter(s) [", paste(paramNames[noSel], collapse = ", "),
                "] unknown and discarded.")

    if (MTD) {
        mtdParams <- formals(makeMTD)
        mtdParams[paramNames[selMTD]] <- params[selMTD]
        formals(makeMTD) <- mtdParams
        mtd <- makeMTD()
    }

    if (what == "PEP") {
        pepParams <- formals(makePEP)
        pepParams[paramNames[selPEP]] <- params[selPEP]
        formals(makePEP) <- pepParams
        tab <- makePEP()
        PRE1 <- "PEH"
        PRE2 <- "PEP"
    } else { ## PRT
        prtParams <- formals(makePRT)
        prtParams[paramNames[selPRT]] <- params[selPRT]
        formals(makePRT) <- prtParams
        tab <- makePRT()
        PRE1 <- "PRH"
        PRE2 <- "PRT"
    }

    if (MTD) 
        cat(mtd, file = file, append = append, sep = "")  
    cat(paste0(PRE1, "\t", paste(colnames(tab), collapse = "\t"), "\n"), sep = "", file = file, append = TRUE)  
    cat(paste0(PRE2, "\t", apply(tab, 1, paste, collapse = "\t"), "\n"), sep = "", file = file, append = TRUE)
    invisible(NULL)
}

##' This function can be used to create a \code{"\linkS4class{MSnSet}"}
##' by reading and parsing an \code{mzTab} file. The metadata section
##' is always used to populate the \code{MSnSet}'s \code{experimentData}
##' slot. 
##'
##' @title Read an 'mzTab' file
##' @param file A \code{character} with the \code{mzTab} file to
##' be read in.
##' @param what One of \code{"PRT"} or \code{"PEP"}, defining
##' which of protein of peptide section should be parse. The metadata
##' section, when available, is always used to populate the
##' \code{experimentData} slot.
##' @param verbose Produce verbose output.
##' @return An instance of class \code{MSnSet}.
##' @author Laurent Gatto
##' @seealso \code{\link{writeMzTabData}} to save an
##' \code{"\linkS4class{MSnSet}"} as an \code{mzTab} file.
##' @examples
##' testfile <- "http://mztab.googlecode.com/svn/legacy/jmztab-1.0/examples/mztab_itraq_example.txt"
##' prot <- readMzTabData(testfile, "PRT")
##' prot
##' pep <- readMzTabData(testfile, "PEP")
##' pep
readMzTabData <- function(file,
                          what = c("PRT", "PEP"),
                          verbose = TRUE) {
    warning("Support for mzTab version 0.9 only. Support will be added soon.")
    
    what <- match.arg(what)
    
    ## .parse1 <- function(x)
    ##   sapply(x, function(.x) strsplit(.x, "\t")[[1]][-1])
    .parse <- function(x) {    
        x <- sapply(x, function(.x) sub("\t", ":", .x))
        names(x) <- sub("^(.+sub\\[[0-9]*\\]).+$", "\\1", x, perl = TRUE)
        x[order(names(x))]
    }
    
    lns <- readLines(file)
    ans <- new("MSnSet")

    ## metadata section
    mtd <- grep("^MTD", lns, value = TRUE) 
    if (length(mtd) > 0) {
        if (verbose)
            message("Detected a metadata section")
        mtd <- sub("MTD\t", "", mtd)
        if (length(title <- grep("-title\t", mtd, value = TRUE)) > 0)
            ans@experimentData@title <- .parse(title)
        if (length(description <- grep("-description\t", mtd, value = TRUE)) > 0) {
            ans@experimentData@other$description <- .parse(description) 
            ## description <- description[-grep("sub", description)] ## added later to experimentData@samples
            ## if (length(description) > 0)
            ##   ans@experimentData@other$description <- .parse(description)
        }
        if (length(sampleProcessing <- grep("-sample_processing\t", mtd, value = TRUE)) > 0)
            ans@experimentData@other$sampleProcessing <- .parse(sampleProcessing)
        if (length(instrumentSource <- grep("-instrument\\[[0-9]*\\]-source\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@ionSource <- .parse(instrumentSource)
        if (length(instrumentAnalyzer <- grep("-instrument\\[[0-9]*\\]-analyzer\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@analyser <- .parse(instrumentAnalyzer)
        if (length(instrumentDetector <- grep("-instrument\\[[0-9]*\\]-detector\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@detectorType <- .parse(instrumentDetector)
        if (length(software <- grep("-software\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@preprocessing <- as.list(.parse(software))
        if (length(fdr <- grep("-false_discovery_rate\t", mtd, value = TRUE)) > 0)
            ans@experimentData@other$fdr <- .parse(fdr)
        if (length(publications <- grep("-publication\t", mtd, value = TRUE)) > 0)
            ans@experimentData@pubMedIds <- .parse(publications)
        if (length(name <- grep("-contact\\[[0-9]*\\]-name\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@name <- .parse(name)
        if (length(affiliation <- grep("-contact\\[[0-9]*\\]-affiliation\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@lab <- .parse(affiliation)
        if (length(email <- grep("-contact\\[[0-9]*\\]-email\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@email <- .parse(email)
        if (length(uri <- grep("-uri\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@url <- .parse(uri)
        if (length(mod <- grep("-mod\t", mtd, value = TRUE)) > 0)
            ans@experimentData@other$modifications <- .parse(mod)
        if (length(modprob <- grep("-mod-probability_method\t", mtd, value = TRUE)) > 0)
            ans@experimentData@other$modProbabilityMethod <- .parse(modprob)
        if (length(quantMethod <- grep("-quantification_method\t", mtd, value = TRUE)) > 0)
            ans@experimentData@other$quantificationMethod <- .parse(quantMethod)
        if (length(protQuantUnit <- grep("-protein-quantification_unit\t", mtd, value = TRUE)) > 0)
            ans@experimentData@other$protQuantUnit <- .parse(protQuantUnit)
        if (length(pepQuantUnit <- grep("-peptide-quantification_unit\t", mtd, value = TRUE)) > 0)
            ans@experimentData@other$pepQuantUnit <- .parse(pepQuantUnit)
        if (length(msFileFormat <- grep("-ms_file\\[[0-9]*\\]-format\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@other$msFileFormat <- .parse(msFileFormat)
        if (length(msFile <- grep("-ms_file\\[[0-9]*\\]-location\t", mtd, value = TRUE)) > 0) 
            ans@processingData@files <- .parse(msFile)
        if (length(msFileIdFormat <- grep("-ms_file\\[[0-9]*\\]-id_format\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@other$msFileIdFormat <- .parse(msFileIdFormat)
        if (length(custom <- grep("-custom\t", mtd, value = TRUE)) > 0) {
            custom <- custom[-grep("sub", custom)] ## added later to experimentData@samples
            if (length(custom) > 0)
                ans@experimentData@other$custom <- .parse(custom)
        }
        ## sub metadata
        if (length(species <- grep("-species\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@samples <- append(ans@experimentData@samples, list(species = .parse(species)))
        if (length(tissue <- grep("-tissue\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0)      
            ans@experimentData@samples <- append(ans@experimentData@samples, list(tissue = .parse(tissue)))
        if (length(cellType <- grep("-cell_type\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0)      
            ans@experimentData@samples <- append(ans@experimentData@samples, list(cellType = .parse(cellType)))
        if (length(disease <- grep("-disease\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0)      
            ans@experimentData@samples <- append(ans@experimentData@samples, list(disease = .parse(disease)))
        if (length(description_ <- grep("-sub\\[[0-9]*\\]-description\t", mtd, value = TRUE)) > 0)      
            ans@experimentData@samples <- append(ans@experimentData@samples, list(description = .parse(description_)))
        if (length(quantReagent_ <- grep("-sub\\[[0-9]*\\]-quantification_reagent\t", mtd, value = TRUE)) > 0)      
            ans@experimentData@samples <- append(ans@experimentData@samples, list(quantReagent = .parse(quantReagent_)))
        if (length(custom_ <- grep("-sub\\[[0-9]*\\]-custom\t", mtd, value = TRUE)) > 0) 
            ans@experimentData@samples <- append(ans@experimentData@samples, list(custom = .parse(custom_)))
    }

    if (what == "PRT") {
        prh <- grep("^PRH", lns, value = TRUE)
        prt <- grep("^PRT", lns, value = TRUE)
        if (length(prh) != 1 & length(prt) < 1)
            stop("No protein section found.")
        if (verbose)
            message("Detected a protein section")
        prt <- sub("^PRT\t", "", prt)
        prh <- sub("^PRH\t", "", prh)
        l <- sapply(prt, strsplit, "\t")
        names(l) <- NULL
        tab <- data.frame(do.call(rbind, l))
        tabnms <- strsplit(prh, "\t")[[1]]
        tabnms <- tabnms[tabnms != ""]
        names(tab) <- tabnms 
    } else { ## "PEP"
        peh <- grep("^PEH", lns, value = TRUE)
        pep <- grep("^PEP", lns, value = TRUE)
        if (length(peh) != 1 & length(pep) < 1)
            stop("No peptide section found.")
        if (verbose)
            message("Detected a peptide section")
        pep <- sub("^PEP\t", "", pep)
        peh <- sub("^PEH\t", "", peh)
        l <- sapply(pep, strsplit, "\t")
        names(l) <- NULL
        tab <- data.frame(do.call(rbind, l))
        tabnms <- strsplit(peh, "\t")[[1]]
        tabnms <- tabnms[tabnms != ""]
        names(tab) <- tabnms 
    }
    
    if (any(esetCols <- grepl("_abundance_sub\\[[0-9]*\\]", names(tab)))) {
        eset <- as.matrix(tab[, esetCols])
        mode(eset) <- "numeric"
        nms <- colnames(eset)
        colnames(eset) <- sub("^.+_abundance_", "", colnames(eset))    
        fdata <- new("AnnotatedDataFrame", data = tab[, !esetCols])
        pdata <- new("AnnotatedDataFrame",
                     data = data.frame(abundance = nms, row.names = colnames(eset)))
    } else {
        warning("No quantitative data found - empty assayData slot.")
        eset <- matrix(nrow = nrow(tab), ncol = 0)
        fdata <- new("AnnotatedDataFrame", data = tab)
        pdata <- new("AnnotatedDataFrame")
    }

    rownames(eset) <- featureNames(fdata)
    exprs(ans) <- eset
    featureData(ans) <- fdata
    phenoData(ans) <- pdata

    ## adding mzTab file
    ans@processingData@files <-
        c(ans@processingData@files,
          file)
    ans@processingData@processing <-
        c(ans@processingData@processing,
          paste0("mzTab read: ", date()))
    
    if (validObject(ans))
        return(ans)
}
