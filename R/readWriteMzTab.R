##' This function can be used to create an
##' \code{"\linkS4class{MSnSet}"} by reading and parsing an
##' \code{mzTab} file. The metadata section is always used to populate
##' the \code{MSnSet}'s \code{experimentData()@@other$mzTab} slot.
##'
##' @title Read an 'mzTab' file
##'
##' @param file A \code{character} with the \code{mzTab} file to be
##'     read in.
##'
##' @param what One of \code{"PRT"}, \code{"PEP"} or \code{"PSM"},
##'     defining which of protein, peptide PSMs section should be
##'     returned as an \code{MSnSet}.
##'
##' @param version A \code{character} defining the format
##'     specification version of the mzTab file. Default is
##'     \code{"1.0"}. Version \code{"0.9"} is available of backwards
##'     compatibility. See \code{\link{readMzTabData_v0.9}} for
##'     details.
##'
##' @param verbose Produce verbose output.
##'
##' @return An instance of class \code{MSnSet}.
##'
##' @seealso See \code{\link{MzTab}} and \code{\link{MSnSetList}} for
##'     details about the inners of \code{readMzTabData}.
##'
##' @author Laurent Gatto
##' @examples
##' testfile <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/1_0-Proteomics-Release/PRIDE_Exp_Complete_Ac_16649.xml-mztab.txt"
##'
##' prot <- readMzTabData(testfile, "PRT")
##'
##' prot
##'
##' head(fData(prot))
##'
##' head(exprs(prot))
##'
##' psms <- readMzTabData(testfile, "PSM")
##'
##' psms
##'
##' head(fData(psms))
readMzTabData <- function(file, what = c("PRT", "PEP", "PSM"),
                          version = c("1.0", "0.9"),
                          verbose = isMSnbaseVerbose()) {
    version <- match.arg(version)
    what <- match.arg(what)
    if (version == "0.9") {
        if (what == "PSM") stop("Only 'PRT' or 'PEP' supported in mzTab version 0.9.")
        readMzTabData_v0.9(file, what, verbose)
    } else {
        ans <- as(MzTab(file), "MSnSetList")
        ans <- switch(what,
                      PRT = ans[["Proteins"]],
                      PEP = ans[["Peptides"]],
                      PSM = ans[["PSMs"]])
        return(ans)
    }
}

##' @title Export an MzTab object as mzTab file.
##'
##' @description
##'
##' `writeMzTabData` exports an [MzTab] object as mzTab file.  Note
##' that the comment section "COM" are not written out.
##'
##' @param object [MzTab] object, either read in by `MzTab()` or assembled.
##'
##' @param file `character(1)` with the file name.
##'
##' @param what `character` with names of the sections to be written
##'     out. Expected sections are `"MT"`, `"PEP"`, `"PRT"`, `"PSM"`,
##'     `"SML"`, `"SMF"`, or `"SME"`.
##'
##' @author Steffen Neumann
##'
##' @md
writeMzTabData <- function(object, file,
                           what = c("MT", "PEP", "PRT", "PSM", "SML", "SMF", "SME")) {
    if ("MT" %in% what) {
        mtdsection <- cbind("MTD", names(object@Metadata), unlist(object@Metadata))
        write.table(mtdsection, file = file, append = FALSE, sep = "\t", na = "null", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
    }

    if ("PRT" %in% what && nrow(object@Proteins) > 0) {
        cat("\n", file = file, append = TRUE)
        section <- data.frame("PRH" = "PRT", object@Proteins, check.names = FALSE)
        ## The suppressed warning is: In write.table(..., append = TRUE, col.names = TRUE, ...)
        ## appending column names to file => That's exactly what we want to do in mzTab!
        suppressWarnings(write.table(section, file = file, append = TRUE,
                                     sep = "\t", na = "null", quote = FALSE,
                                     row.names = FALSE, col.names = TRUE))
    }

    if ("PEP" %in% what && nrow(object@Peptides) > 0) {
        cat("\n", file = file, append = TRUE)
        section <- data.frame("PEH" = "PEP", object@Peptides, check.names = FALSE)
        ## The suppressed warning is: In write.table(..., append = TRUE, col.names = TRUE, ...)
        ## appending column names to file => That's exactly what we want to do in mzTab!
        suppressWarnings(write.table(section, file = file, append = TRUE,
                                     sep = "\t", na = "null", quote = FALSE,
                                     row.names = FALSE, col.names = TRUE))
    }

    if ("PSM" %in% what && nrow(object@PSMs) > 0) {
        cat("\n", file = file, append = TRUE)
        section <- data.frame("PSH" = "PSM", object@PSMs, check.names = FALSE)
        ## The suppressed warning is: In write.table(..., append = TRUE, col.names = TRUE, ...)
        ## appending column names to file => That's exactly what we want to do in mzTab!
        suppressWarnings(write.table(section, file = file, append = TRUE,
                                     sep = "\t", na = "null", quote = FALSE,
                                     row.names = FALSE, col.names = TRUE))
    }

    if ("SML" %in% what && nrow(object@SmallMolecules) > 0) {
        cat("\n", file = file, append = TRUE)
        section <- data.frame("SMH" = "SML", object@SmallMolecules, check.names = FALSE)
        ## The suppressed warning is: In write.table(..., append = TRUE, col.names = TRUE, ...)
        ## appending column names to file => That's exactly what we want to do in mzTab!
        suppressWarnings(write.table(section, file = file, append = TRUE,
                                     sep = "\t", na = "null", quote = FALSE,
                                     row.names = FALSE, col.names = TRUE))
    }

    if ("SMF" %in% what && nrow(object@MoleculeFeatures) > 0) {
        cat("\n", file = file, append = TRUE)
        section <- data.frame("SFH" = "SMF", object@MoleculeFeatures, check.names = FALSE)
        ## The suppressed warning is: In write.table(..., append = TRUE, col.names = TRUE, ...)
        ## appending column names to file => That's exactly what we want to do in mzTab!
        suppressWarnings(write.table(section, file = file, append = TRUE,
                                     sep = "\t", na = "null", quote = FALSE,
                                     row.names = FALSE, col.names = TRUE))
    }

    if ("SME" %in% what && nrow(object@MoleculeEvidence) > 0) {
        cat("\n", file = file, append = TRUE)
        section <- data.frame("SEH" = "SME", object@MoleculeEvidence, check.names = FALSE)
        ## The suppressed warning is: In write.table(..., append = TRUE, col.names = TRUE, ...)
        ## appending column names to file => That's exactly what we want to do in mzTab!
        suppressWarnings(write.table(section, file = file, append = TRUE,
                                     sep = "\t", na = "null", quote = FALSE,
                                     row.names = FALSE, col.names = TRUE))
    }
}

## ===================================
## Legacy code

makeMTD <- function(...) .Defunct()
makePEP <- function(...) .Defunct()
makePRT <- function(...) .Defunct()

##' This function can be used to create a \code{"\linkS4class{MSnSet}"}
##' by reading and parsing an \code{mzTab} file. The metadata section
##' is always used to populate the \code{MSnSet}'s \code{experimentData}
##' slot.
##'
##' @title Read an 'mzTab' file
##'
##' @param file A \code{character} with the \code{mzTab} file to be
##'     read in.
##'
##' @param what One of \code{"PRT"} or \code{"PEP"}, defining which of
##'     protein of peptide section should be parse. The metadata
##'     section, when available, is always used to populate the
##'     \code{experimentData} slot.
##'
##' @param verbose Produce verbose output.
##'
##' @return An instance of class \code{MSnSet}.
##'
##' @author Laurent Gatto
##'
##' @seealso \code{\link{writeMzTabData}} to save an
##'     \code{"\linkS4class{MSnSet}"} as an \code{mzTab} file.
##'
readMzTabData_v0.9 <- function(file,
                               what = c("PRT", "PEP"),
                               verbose = isMSnbaseVerbose())
    .Defunct(msg = "Version 0.9 is deprecated. Please see '?readMzTabData' and '?MzTab' for details.")
