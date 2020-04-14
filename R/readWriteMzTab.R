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
##' @examples
##' testfile <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/legacy/jmztab-1.0/examples/mztab_itraq_example.txt"
##' 
##' prot <- readMzTabData_v0.9(testfile, "PRT")
##' 
##' prot
##' 
##' pep <- readMzTabData_v0.9(testfile, "PEP")
##' 
##' pep
readMzTabData_v0.9 <- function(file,
                               what = c("PRT", "PEP"),
                               verbose = isMSnbaseVerbose()) {
    .Deprecated(msg = "Version 0.9 is deprecated. Please see '?readMzTabData' and '?MzTab' for details.")
    what <- match.arg(what)
    
    ## .parse1 <- function(x)
    ##   sapply(x, function(.x) strsplit(.x, "\t")[[1]][-1])
    .parse <- function(x) {
        x <- sapply(x, function(.x) sub("\t", ":", .x))
        names(x) <- sub("^(.+sub\\[[0-9]*\\]).+$", "\\1", x, perl = TRUE)
        x[order(names(x))]
    }

    processingDataFiles <- NULL
    
    lns <- readLines(file)
    miape <- new("MIAPE")

    ## metadata section
    mtd <- grep("^MTD", lns, value = TRUE)
    if (length(mtd) > 0) {
        if (verbose)
            message("Detected a metadata section")
        mtd <- sub("MTD\t", "", mtd)
        if (length(title <- grep("-title\t", mtd, value = TRUE)) > 0)
            miape@title <- .parse(title)
        if (length(description <- grep("-description\t", mtd, value = TRUE)) > 0) {
            miape@other$description <- .parse(description) 
            ## description <- description[-grep("sub", description)] ## added later to experimentData@samples
            ## if (length(description) > 0)
            ##   miape@other$description <- .parse(description)
        }
        if (length(sampleProcessing <- grep("-sample_processing\t", mtd, value = TRUE)) > 0)
            miape@other$sampleProcessing <- .parse(sampleProcessing)
        if (length(instrumentSource <- grep("-instrument\\[[0-9]*\\]-source\t", mtd, value = TRUE)) > 0) 
            miape@ionSource <- .parse(instrumentSource)
        if (length(instrumentAnalyzer <- grep("-instrument\\[[0-9]*\\]-analyzer\t", mtd, value = TRUE)) > 0) 
            miape@analyser <- .parse(instrumentAnalyzer)
        if (length(instrumentDetector <- grep("-instrument\\[[0-9]*\\]-detector\t", mtd, value = TRUE)) > 0) 
            miape@detectorType <- .parse(instrumentDetector)
        if (length(software <- grep("-software\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0) 
            miape@preprocessing <- as.list(.parse(software))
        if (length(fdr <- grep("-false_discovery_rate\t", mtd, value = TRUE)) > 0)
            miape@other$fdr <- .parse(fdr)
        if (length(publications <- grep("-publication\t", mtd, value = TRUE)) > 0)
            miape@pubMedIds <- .parse(publications)
        if (length(name <- grep("-contact\\[[0-9]*\\]-name\t", mtd, value = TRUE)) > 0) 
            miape@name <- .parse(name)
        if (length(affiliation <- grep("-contact\\[[0-9]*\\]-affiliation\t", mtd, value = TRUE)) > 0) 
            miape@lab <- .parse(affiliation)
        if (length(email <- grep("-contact\\[[0-9]*\\]-email\t", mtd, value = TRUE)) > 0) 
            miape@email <- .parse(email)
        if (length(uri <- grep("-uri\t", mtd, value = TRUE)) > 0) 
            miape@url <- .parse(uri)
        if (length(mod <- grep("-mod\t", mtd, value = TRUE)) > 0)
            miape@other$modifications <- .parse(mod)
        if (length(modprob <- grep("-mod-probability_method\t", mtd, value = TRUE)) > 0)
            miape@other$modProbabilityMethod <- .parse(modprob)
        if (length(quantMethod <- grep("-quantification_method\t", mtd, value = TRUE)) > 0)
            miape@other$quantificationMethod <- .parse(quantMethod)
        if (length(protQuantUnit <- grep("-protein-quantification_unit\t", mtd, value = TRUE)) > 0)
            miape@other$protQuantUnit <- .parse(protQuantUnit)
        if (length(pepQuantUnit <- grep("-peptide-quantification_unit\t", mtd, value = TRUE)) > 0)
            miape@other$pepQuantUnit <- .parse(pepQuantUnit)
        if (length(msFileFormat <- grep("-ms_file\\[[0-9]*\\]-format\t", mtd, value = TRUE)) > 0) 
            miape@other$msFileFormat <- .parse(msFileFormat)
        if (length(msFile <- grep("-ms_file\\[[0-9]*\\]-location\t", mtd, value = TRUE)) > 0) 
            processingDataFiles <- .parse(msFile)
        if (length(msFileIdFormat <- grep("-ms_file\\[[0-9]*\\]-id_format\t", mtd, value = TRUE)) > 0) 
            miape@other$msFileIdFormat <- .parse(msFileIdFormat)
        if (length(custom <- grep("-custom\t", mtd, value = TRUE)) > 0) {
            custom <- custom[-grep("sub", custom)] ## added later to experimentData@samples
            if (length(custom) > 0)
                miape@other$custom <- .parse(custom)
        }
        ## sub metadata
        if (length(species <- grep("-species\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0) 
            miape@samples <- append(miape@samples, list(species = .parse(species)))
        if (length(tissue <- grep("-tissue\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0)      
            miape@samples <- append(miape@samples, list(tissue = .parse(tissue)))
        if (length(cellType <- grep("-cell_type\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0)      
            miape@samples <- append(miape@samples, list(cellType = .parse(cellType)))
        if (length(disease <- grep("-disease\\[[0-9]*\\]\t", mtd, value = TRUE)) > 0)      
            miape@samples <- append(miape@samples, list(disease = .parse(disease)))
        if (length(description_ <- grep("-sub\\[[0-9]*\\]-description\t", mtd, value = TRUE)) > 0)      
            miape@samples <- append(miape@samples, list(description = .parse(description_)))
        if (length(quantReagent_ <- grep("-sub\\[[0-9]*\\]-quantification_reagent\t", mtd, value = TRUE)) > 0)      
            miape@samples <- append(miape@samples, list(quantReagent = .parse(quantReagent_)))
        if (length(custom_ <- grep("-sub\\[[0-9]*\\]-custom\t", mtd, value = TRUE)) > 0) 
            miape@samples <- append(miape@samples, list(custom = .parse(custom_)))
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
        eset[eset == "null"] <- NA
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
    colnames(eset) <- featureNames(pdata)

    ans <- new("MSnSet",
               exprs = eset,
               featureData = fdata,
               phenoData = pdata,
               experimentData = miape)

    ## adding mzTab file
    ans@processingData@files <-
        c(ans@processingData@files,
          processingDataFiles, ## possibly NULL
          file)
    ans@processingData@processing <-
        c(ans@processingData@processing,
          paste0("mzTab read: ", date()))
    
    if (validObject(ans))
        return(ans)
}
