.MzTabList <- setClass("MzTabList",
                       slots = c(
                           Metadata = "list",
                           Proteins = "data.frame",
                           Peptides = "data.frame",
                           Spectra = "data.frame",
                           SmallMolecules = "data.frame",
                           Comments = "character"))

setMethod("show", "MzTabList",
          function(object) {
              cat("Object of class \"", class(object),"\".\n", sep = "")
              descr <- paste0(" Description: ", object@Metadata$description, "\n")
              descr <- paste0(" ", strwrap(descr), "\n")
              cat(descr, sep = "")              
              cat(" Mode:", object@Metadata$`mzTab-mode`, "\n")
              cat(" Type:", object@Metadata$`mzTab-type`, "\n")
              cat(" Available data: ")
              avbl <- sapply(slotNames(mz)[2:5], function(x) nrow(slot(mz, x)) > 0)
              cat(paste(names(avbl)[which(avbl)], collape = ""), "\n")
          })

## Accessors
setMethod("metadata", "MzTabList",
          function(x, ...) x@Metadata)

setMethod("proteins", "MzTabList",
          function(object, ...) object@Proteins)

setMethod("peptides", "MzTabList",
          function(object, ...) object@Peptides)

setMethod("spectra", "MzTabList",
          function(object, ...) object@Spectra)

smallMolecules <- function(x) x@SmallMolecules

## Notes
## Mandatory header fields
## mzTab-version
## mzTab-type: Identification Quantitation
## mzTab-mode: Summary Complete
## description
## ms_run-location[1-n]
## also
##
## “protein_search_engine_score[1-n]”,
## ”peptide_search_engine_score[1-n]”, “psm_search_engine_score[1-n]”
## and “smallmolecule_search_engine_score[1-n]” MUST be reported for
## every search engine score reported in the corresponding section.
##
## “fixed_mod[1-n]” and “variable_mod [1-n]” MUST be reported. If no
## modifications were searched, specific CV parameters need to be used
## (see Section 5.8).

## Niels:
##
## - The file you sent me is not a valif mzTab file

## library(data.table)
## library(stringi)

MzTabList <- function(file) {
    ## Like readLines, but way faster
    ## lines <- stri_read_lines(file)
    lines <- readLines(file)
    lines <- lines[nzchar(lines)]

    ## Split on the first two characters (so headers stay in
    ## the same group as table content rows)
    lineType <- substring(lines, 1, 2)

    ## Could be stricter in the type checking to check that all
    ## three of the first characters match the 10 allowed types
    ## but since it doesn't affect parsing, I don't think it's
    ## worth bothering.
    allowed_types <- c("CO", "MT", "PR", "PE", "PS", "SM")
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
            linesByType[c("MT", "PR", "PE", "PS", "SM")],
            function(x) {
                if (length(x) == 0) return(data.frame())
                return(read.delim(text = x,
                                  na.strings = c("", "null"),
                                  stringsAsFactors = FALSE)[,-1])
            }),
        c("Metadata", "Proteins", "Peptides",
          "Spectra", "SmallMolecules"))
    
    res[["Metadata"]] <- reshapeMetadata(res[["Metadata"]])

    .MzTabList(Metadata = res[["Metadata"]],
               Proteins = res[["Proteins"]],
               Peptides = res[["Peptides"]],
               Spectra = res[["Spectra"]],
               SmallMolecules = res[["SmallMolecules"]],
               Comments = comments)

}

##' @param x A character vector.
##' @return A data table.
##' @importFrom data.table fread setDT
readTab <- function(x) {
    ## Passing data to fread as a string doesn't work for < 2 lines, so
    ## we have to handle these cases separately.
    if (length(x) == 0) {
        return(data.table())
    }
    if (length(x) == 1) {
        ## [[-1]] is because we no longer care about the identifier
        ## at the start of the line.
        return(setDT(read.delim(text = x)[, -1]))
    }
    ## It isn't clear from the documentation which values are allowed
    ## to represent missing values.  There were cases of "" and "null"
    ## in the examples, but this may need to become an argument to
    ## readMzTabData.
    ## Also, some fiddling with converting columns types may be needed.
    ## See https://github.com/Rdatatable/data.table/issues/729
    fread(paste0(x, collapse = "\n"), sep = "\t", na.strings = c("", "null"))[, -1, with = FALSE]
}

##' @param mtd A \code{data.frame} with 2 columns
##' @return A named list, where each element is either a
##' string or a named character vector of length 4.
reshapeMetadata <- function(mtd) {
    stopifnot(ncol(mtd) >= 2)
    ## Second column of mtd contains keys
    ## Third column contains values
    ## Store results in a list, with keys as names
    metadata <- setNames(vector("list", nrow(mtd)), mtd[[1]])
    metadata[1:length(metadata)] <- mtd[[2]]
    
    ## ## Need different behaviour for strings vs. controlled vocab
    ## rx <- "\\[([[:alnum:]]+), *([[:alnum:]]+:[[:digit:]]+), *([[:print:]]+), *([[:print:]]+)?\\]"
    ## isControlledVocab <- stri_detect_regex(mtd[[2]], rx)
    ## controlledVocabValues <- stri_match_all_regex(
    ##     mtd[[2]][isControlledVocab],
    ##     rx
    ##     )
    ## metadata[!isControlledVocab] <- mtd[[2]][!isControlledVocab]
    ## metadata[isControlledVocab] <- lapply(
    ##     controlledVocabValues,
    ##     function(x)
    ##         {
    ##             setNames(x[-1], c("vocab", "id", "name", "def"))
    ##         }
    ##     )
    metadata
}


## testing
allmzt <- dir("~/dev/00_github/MSnbase/sandbox/mzTabExamples/", full.names=TRUE)
sapply(allmzt, MzTabList)

## coerce MzTabList as MSnSetList|MSnset

## what about readMzTabData? Keep it for backwards compatibility?
## deprecate writeMzTabData
