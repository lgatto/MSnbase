setMethod("show", "MzTab",
          function(object) {
              cat("Object of class \"", class(object),"\".\n", sep = "")
              descr <- paste0(" Description: ", object@Metadata$description, "\n")
              descr <- paste0(" ", strwrap(descr), "\n")
              cat(descr, sep = "")
              cat(" Mode:", object@Metadata$`mzTab-mode`, "\n")
              cat(" Type:", object@Metadata$`mzTab-type`, "\n")
              cat(" Available data: ")
              avbl <- sapply(slotNames(object)[3:6],
                             function(x) nrow(slot(object, x)) > 0)
              cat(paste(names(avbl)[which(avbl)], collape = ""), "\n")
          })

#' @title Export an MzTab object as mzTab file.
#'
#' @description
#'
#' `write.mzTab` exports an MzTab object as mzTab file. 
#' Note that the comment section "COM" are not written out. 
#'
#' @param object [MzTab] object, either read in by MzTab() or assembled. 
#'
#' @param file `character` with the file name. 
#'
#' @param what `character` with names of the sections to be written out. 
#'
#' @author Steffen Neumann
#'
#' @md
#'
setMethod("write.mzTab", 
          "MzTab",
          function(object,
                   file, 
                   what = c("MT", "PEP", "PRT", "PSM", "SML", "SMF", "SME"),
                   ...) {
            
            if ("MT" %in% what) {
              mtdsection <- cbind("MTD", names(object@Metadata), unlist(object@Metadata))
              write.table(mtdsection, file=file, append=FALSE, sep="\t", na="null", quote=FALSE, 
                          row.names=FALSE, col.names=FALSE)
            }
            
            if ("PRT" %in% what && nrow(object@Proteins) >0) {
              cat("\n", file=file, append=TRUE)    
              section <- data.frame("PRH"="PRT", object@Proteins, check.names = FALSE)
              ## The supporessed warning is: In write.table(..., append = TRUE, col.names=TRUE, ...)
              ## appending column names to file => That's exactly what we want to do in mzTab!
              suppressWarnings(write.table(section, file=file, append=TRUE, sep="\t", na="null", quote=FALSE, 
                                           row.names=FALSE, col.names=TRUE))
            }
            
            if ("PEP" %in% what && nrow(object@Peptides) >0) {
              cat("\n", file=file, append=TRUE)    
              section <- data.frame("PEH"="PEP", object@Peptides, check.names = FALSE)
              ## The supporessed warning is: In write.table(..., append = TRUE, col.names=TRUE, ...)
              ## appending column names to file => That's exactly what we want to do in mzTab!
              suppressWarnings(write.table(section, file=file, append=TRUE, sep="\t", na="null", quote=FALSE, 
                                           row.names=FALSE, col.names=TRUE))
            }
            
            if ("PSM" %in% what && nrow(object@PSMs) >0) {
              cat("\n", file=file, append=TRUE)    
              section <- data.frame("PSH"="PSM", object@PSMs, check.names = FALSE)
              ## The supporessed warning is: In write.table(..., append = TRUE, col.names=TRUE, ...)
              ## appending column names to file => That's exactly what we want to do in mzTab!
              suppressWarnings(write.table(section, file=file, append=TRUE, sep="\t", na="null", quote=FALSE, 
                                           row.names=FALSE, col.names=TRUE))
            }
            
            if ("SML" %in% what && nrow(object@SmallMolecules) >0) {
              cat("\n", file=file, append=TRUE)    
              section <- data.frame("SMH"="SML", object@SmallMolecules, check.names = FALSE)
              ## The supporessed warning is: In write.table(..., append = TRUE, col.names=TRUE, ...)
              ## appending column names to file => That's exactly what we want to do in mzTab!
              suppressWarnings(write.table(section, file=file, append=TRUE, sep="\t", na="null", quote=FALSE, 
                                           row.names=FALSE, col.names=TRUE))
            }
            
            if ("SMF" %in% what && nrow(object@MoleculeFeatures) >0) {
              cat("\n", file=file, append=TRUE)    
              section <- data.frame("SFH"="SMF", object@MoleculeFeatures, check.names = FALSE)
              ## The supporessed warning is: In write.table(..., append = TRUE, col.names=TRUE, ...)
              ## appending column names to file => That's exactly what we want to do in mzTab!
              suppressWarnings(write.table(section, file=file, append=TRUE, sep="\t", na="null", quote=FALSE, 
                                           row.names=FALSE, col.names=TRUE))
            }
            
            if ("SME" %in% what && nrow(object@MoleculeEvidence) >0) {
              cat("\n", file=file, append=TRUE)    
              section <- data.frame("SEH"="SME", object@MoleculeEvidence, check.names = FALSE)
              ## The supporessed warning is: In write.table(..., append = TRUE, col.names=TRUE, ...)
              ## appending column names to file => That's exactly what we want to do in mzTab!
              suppressWarnings(write.table(section, file=file, append=TRUE, sep="\t", na="null", quote=FALSE, 
                                           row.names=FALSE, col.names=TRUE))
            }
          })

## Accessors

## Generic from S4Vectors
setMethod("metadata", "MzTab",
          function(x, ...) x@Metadata)

mzTabMode <- function(x)
    metadata(x)$`mzTab-mode`

mzTabType <- function(x)
    metadata(x)$`mzTab-type`

## Generic from BiocGenerics
setMethod("fileName", "MzTab",
          function(object, ...) object@Filename)

setMethod("fileNames", "MzTab",
          function(object, ...) object@Filename)

## Generic from ProtGenerics
setMethod("proteins", "MzTab",
          function(object, ...) object@Proteins)

## Generic from ProtGenerics
setMethod("peptides", "MzTab",
          function(object, ...) object@Peptides)

## Generic from ProtGenerics
setMethod("psms", "MzTab",
          function(object, ...) object@PSMs)

smallMolecules <- function(x) x@SmallMolecules
moleculeFeatures <- function(x) x@MoleculeFeatures
moleculeEvidence <- function(x) x@MoleculeEvidence

comments <- function(x) x@Comments

## Constructor
##  Based on @richierocks contribution
##  https://github.com/lgatto/MSnbase/issues/41
MzTab <- function(file) {
    file <- file[1]
    lines <- readLines(file)
    lines <- lines[-grep("^\\s*$", readLines(file))]
    lines <- lines[nzchar(lines)]

    ## Split on the first characters, make sure headers stay in
    ## the same group as table content rows
    
    lineType <- sapply(substring(lines, 1, 3), function(s) switch(s,
                       MTD = "MT",
                       COM = "CO",
                       PRH = "PR",
                       PRT = "PR",
                       PEH = "PE",
                       PEP = "PE",
                       PSH = "PS",
                       PSM = "PS",
                       SMH = "SM",
                       SML = "SM",
                       SFH = "SF",
                       SMF = "SF",
                       SEH = "SE",
                       SME = "SE"))

    ## Could be stricter in the type checking to check that all
    ## three of the first characters match the 10 allowed types
    ## but since it doesn't affect parsing, I don't think it's
    ## worth bothering.
    allowed_types <- c("CO", "MT", "PR", "PE", "PS", "SM", "SF", "SE")
    stopifnot(all(lineType %in% allowed_types))
    linesByType <- split(lines, lineType)

    ## Comments are easy: just strip the first four characters
   ## from each line.  Though is it important to record the
    ## position in the file where they were found?
    comments <- substring(linesByType[["CO"]], 5)

    ## Parse the other five blocks in a loop, then fix up
    ## metadata afterwards
    res <- setNames(
        lapply(
            linesByType[c("MT", "PR", "PE", "PS", "SM", "SF", "SE")],
            function(x) {
                if (length(x) == 0) return(data.frame())
                return(read.delim(text = x,
                                  header=ifelse(all(grepl("^MTD", x)), FALSE, TRUE), ## MTD has no header
                                  na.strings = c("", "null"),
                                  check.names = FALSE,
                                  stringsAsFactors = FALSE)[,-1])
            }),
        c("Metadata", "Proteins", "Peptides", "PSMs", "SmallMolecules",
          "MoleculeFeatures", "MoleculeEvidence"))
    
    res[["Metadata"]] <- reshapeMetadata(res[["Metadata"]])

    .MzTab(Metadata = res[["Metadata"]],
           Filename = file,
           Proteins = res[["Proteins"]],
           Peptides = res[["Peptides"]],
           PSMs = res[["PSMs"]],
           SmallMolecules = res[["SmallMolecules"]],
           MoleculeFeatures = res[["MoleculeFeatures"]],
           MoleculeEvidence = res[["MoleculeEvidence"]],
           Comments = comments)

}

##' @param mtd A \code{data.frame} with 2 columns
##' @return A named list, where each element is a string
##' @noRd
reshapeMetadata <- function(mtd) {
    stopifnot(ncol(mtd) >= 2)
    metadata <- setNames(vector("list", nrow(mtd)), mtd[[1]])
    metadata[1:length(metadata)] <- mtd[[2]]    
    metadata
}

setAs("MzTab", "MSnSetList",
      function(from, to = "MSnSetList") {
          MSnSetList(list(
              Proteins = makeProtMSnSet(from),
              Peptides = makePepMSnSet(from),
              PSMs = makePsmMSnSet(from)))
      })


##' @param x MSnSet object
##' @param y MzTab object
##' @return x decorated with metadata(y)
##' @noRd
addMzTabMetadata <- function(x, y) {
    experimentData(x)@other$mzTab <- metadata(y)
    if (any(i <- grepl("publication", names(metadata(y))))) 
        pubMedIds(x) <- unlist(metadata(y)[i], use.names = FALSE)
    x@processingData@files <- fileName(y)
    ## This would need www access, if to use rols
    ## experimentData(x)@samples$species <- 'sample[1-n]-species[1-n]'
    if (validObject(x)) x
}

makeProtMSnSet <- function(object,
                           protabundance = "protein_abundance_assay") {
    x <- proteins(object)
    if (nrow(x) == 0) {
        ans <- new("MSnSet")
    } else {
        ecols <- grep(protabundance, names(x))        
        e <- as.matrix(x[, ecols])
        if (length(ecols) > 0) fd <- x[, -ecols]
        else fd <- x
        rownames(e) <- rownames(fd) <-
            make.names(fd[, "accession"], unique = TRUE)
        pd <- data.frame(row.names = colnames(e))
        ans <- MSnSet(exprs = e, fData = fd, pData = pd)
    }
    addMzTabMetadata(ans, object)
}

makePepMSnSet <- function(object,
                          pepabundance = "peptide_abundance_assay") {
    x <- peptides(object)
    if (nrow(x) == 0) {
        ans <- new("MSnSet")
    } else {
        ecols <- grep(pepabundance, names(x))
        e <- as.matrix(x[, ecols])
        if (length(ecols) > 0) fd <- x[, -ecols]
        else fd <- x
        rownames(e) <- rownames(fd) <-
            make.names(fd[, "sequence"], unique = TRUE)
        pd <- data.frame(row.names = colnames(e))
        ans <- MSnSet(exprs = e, fData = fd, pData = pd)
    }
    addMzTabMetadata(ans, object)
}

makePsmMSnSet <- function(object) {
    x <- psms(object)
    if (nrow(x) == 0) {
        ans <- new("MSnSet")
    } else {
        e <- matrix(NA_real_, nrow = nrow(x), ncol = 0)
        fd <- x
        rownames(e) <- rownames(fd) <-
            make.names(fd[, "PSM_ID"], unique = TRUE)
        pd <- data.frame(row.names = colnames(e))
        ans <- MSnSet(exprs = e, fData = fd, pData = pd)        
    }
    addMzTabMetadata(ans, object)
}
