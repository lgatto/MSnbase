## What about a validity method for MzTab objects?
##
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

setMethod("show", "MzTab",
          function(object) {
              cat("Object of class \"", class(object),"\".\n", sep = "")
              descr <- paste0(" Description: ", object@Metadata$description, "\n")
              descr <- paste0(" ", strwrap(descr), "\n")
              cat(descr, sep = "")              
              cat(" Mode:", object@Metadata$`mzTab-mode`, "\n")
              cat(" Type:", object@Metadata$`mzTab-type`, "\n")
              cat(" Available data: ")
              avbl <- sapply(slotNames(object)[2:5],
                             function(x) nrow(slot(object, x)) > 0)
              cat(paste(names(avbl)[which(avbl)], collape = ""), "\n")
          })

## Accessors
## Generic from S4Vectors
setMethod("metadata", "MzTab",
          function(x, ...) x@Metadata)

## Generic from ProtGenerics
setMethod("proteins", "MzTab",
          function(object, ...) object@Proteins)

## Generic from ProtGenerics
setMethod("peptides", "MzTab",
          function(object, ...) object@Peptides)

## Generic from ProtGenerics
setMethod("spectra", "MzTab",
          function(object, ...) object@Spectra)

smallMolecules <- function(x) x@SmallMolecules

comments <- function(x) x@Comments

## Constructor
##  Based on @richierocks contribution
##  https://github.com/lgatto/MSnbase/issues/41
MzTab <- function(file) {
    file <- file[1]
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
        c("Metadata", "Proteins", "Peptides", "Spectra", "SmallMolecules"))
    
    res[["Metadata"]] <- reshapeMetadata(res[["Metadata"]])

    .MzTab(Metadata = res[["Metadata"]],
               Proteins = res[["Proteins"]],
               Peptides = res[["Peptides"]],
               Spectra = res[["Spectra"]],
               SmallMolecules = res[["SmallMolecules"]],
               Comments = comments)

}

##' @param mtd A \code{data.frame} with 2 columns
##' @return A named list, where each element is a string
reshapeMetadata <- function(mtd) {
    stopifnot(ncol(mtd) >= 2)
    metadata <- setNames(vector("list", nrow(mtd)), mtd[[1]])
    metadata[1:length(metadata)] <- mtd[[2]]    
    metadata
}

## coerce MzTab as MSnSetList (peps, prots, spectra) or MSnset (only
## one)

## what about readMzTabData? Keep it for backwards compatibility?
## deprecate writeMzTabData
