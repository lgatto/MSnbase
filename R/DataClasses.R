## Returns of class version as documented in .MSnBaseEnd$ClassVersions
## as and instance of class Versions.
getClassVersion <- function(x) {
    if (!is.character(x))
        x <- class(x)[1]
    ## This get class versions from parent classes (if any)
    ver <- classVersion(x)
    ## Adds (or overwrites) x's class version to the list of class
    ## versions
    ver[x] <- getClassVersionString(x)
    ver
}
## Utility to just extract the version string from the environment.
getClassVersionString <- function(x) {
    if (!is.character(x))
        x <- class(x)[1]
    return(.MSnbaseEnv$ClassVersions[x])
}

######################################################################
## MSnProcess: Container for MSnExp and MSnSet processing information
## See online documentation for more information.
setClass("MSnProcess",
         representation = representation(
             files = "character",
             processing = "character",
             merged = "logical",
             cleaned = "logical",
             removedPeaks = "character",
             smoothed = "logical",
             trimmed = "numeric",
             normalised = "logical",
             MSnbaseVersion = "character"),
         contains = c("Versioned"),
         prototype  =  prototype(
             new("Versioned", versions = c(MSnProcess = "0.1.3")),
             processing = character(),
             files = character(),
             trimmed = numeric(),
             removedPeaks = character(),
             MSnbaseVersion = character())) ## set in initialize()

#################################################################
## The 'Minimum Information About a Proteomics Experiment' Class
## See online documentation for more information.
setClass("MIAPE",
         representation = representation(
             title = "character",
             url = "character",
             ##  Publication details
             abstract = "character",
             pubMedIds = "character",
             ## Other useful slots, from MIAME
             samples = "list",
             preprocessing = "list",
             other = "list",
             ## ########################
             ## Based on MIAPE-MS 2.24
             ##  will be updated with MIAPE-MSI and MIAPE-Quant
             ## 1. General features - (a) Global descriptors
             dateStamp = "character",
             ##    Responsible person
             name = "character",
             lab = "character",
             contact = "character",
             email = "character",
             ##    Instrument details
             instrumentModel = "character",
             instrumentManufacturer = "character",
             instrumentCustomisations = "character",
             ## 1. General features - (b) Control and analysis software
             softwareName = "character",
             softwareVersion = "character",
             switchingCriteria = "character",
             isolationWidth = "numeric",
             parameterFile = "character",
             ## 2. Ion sources -- will be updated to
             ##                   provided details specific to
             ##                   different sources
             ionSource = "character", ## ESI, MALDI, ...
             ionSourceDetails = "character",
             ## 3. Post-source component
             analyser = "character", ## Quad, TOF, Trap, ...
             analyserDetails = "character",
             ## 3. Post-source component - (d) Collision cell
             collisionGas = "character",
             collisionPressure = "numeric",
             collisionEnergy = "character",
             ## 3. Post-source component - (f) Detectors
             detectorType = "character",
             detectorSensitivity = "character"
             ## 4. Spectrum and peak list generation and annotation
             ##    (a) Spectrum description
             ##    (b) Peak list generation
             ##    (c) Quantitation for selected ions
             ## -- see Spectrum and MSnProcess objects
             ),
         contains = c("MIAxE"),
         prototype = prototype(
             new("Versioned",
                 versions = c(classVersion("MIAxE"), MIAPE = "0.2.2")),
             name = "",
             lab = "",
             contact = "",
             title = "",
             abstract = "",
             url = "",
             pubMedIds = "",
             email = "",
             samples = list(),
             preprocessing = list(),
             other = list()))

############################################################################
## NAnnotatedDataFrame: As Biobase's AnnotatedDataFrame, it is composed of
## a data.frame, with annotations about columns named
## in the data slot contained in the metadata slot.
## In addition, it contains a multiplex slot to make explicite that
## the AnnotatedDataFrame is applied to a set of mulitplexed tags.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setClass("NAnnotatedDataFrame",
         representation(multiplex = "numeric",
                        multiLabels = "character"),
         contains = c("AnnotatedDataFrame"),
         prototype = prototype(
              new("Versioned", versions = list(NAnnotatedDataFrame="0.0.3")),
             multiplex = 1,
             multiLabels = "Single run"),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (length(object@multiLabels) != object@multiplex)
                 msg <- validMsg(msg, "Number of multiplex does not match it's labels.")
             if (is.null(msg)) TRUE
             else msg
         })


#############################################################################
## pSet: similarly to eSet but with a focus toward proteomics experiments,
## pSet is a VIRTUAL class containing assay data (typically, one or many
## different sets of spectra), phenotypic data (describing the samples involved
## in the experiment), experimental data (describing the methods and
## protocols used) and feature data (describing the features in the assay).
##
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setClass("pSet",
         representation(assayData = "environment", ## locked environment
                        phenoData = "NAnnotatedDataFrame",
                                        # Filenames,  Fractions, replicates
                                        # file1.mzML, 1        , 1          (n1 spectra)
                                        # file2.mzML, 2        , 1          (n2 spectra)
                                        # file3.mzML, 1        , 2          (n3 spectra)
                                        # file4.mzML, 1        , 2          (n4 spectra)
                                        # ...
                                        # How to link individual spectra in assayData env to phenoData?
                                        # Use the spectrum@fromFile slot
                                        # This phenoData will be expanded to an NAnnotatedDataFrame
                                        # when the MSnSet is created from MSnExp
                        featureData = "AnnotatedDataFrame",
                                        # How to link individual spectra in assayData env to featureData?
                                        # The spectra and the rowNames of featureData will be named X1..(n1+n2+...)
                        experimentData = "MIAxE",
                        protocolData   = "AnnotatedDataFrame",
                        processingData = "MSnProcess",
                        ##onDisk         = "logical",
                        .cache = "environment", ## locked
                        "VIRTUAL"),
         contains = "Versioned",
         prototype = prototype(
             new("Versioned", versions = c(pSet = "0.1.1")),
             assayData = new.env(parent=emptyenv()),
             experimentData = new("MIAPE"),
             phenoData = new("NAnnotatedDataFrame",
                 dimLabels=c("sampleNames", "fileNumbers")),
             featureData = new("AnnotatedDataFrame",
                 dimLabels=c("featureNames", "featureColumns")),
             protocolData = new("AnnotatedDataFrame",
                                dimLabels=c("sampleNames", "sampleColumns"))
             ##,onDisk = FALSE
         )
         )

##################################################################
## Container for MSn Experiments Data and Meta-Data
## See online documentation for more information.
setClass("MSnExp",
         contains=c("pSet"),
         prototype = prototype(
             new("VersionedBiobase",
                 versions = c(classVersion("pSet"), MSnExp="0.3.1")),
             experimentData = new("MIAPE"))
         )

###################################################################
## The Spectrum class and it's sub-classes Spectrum1 and Spectrum2
## See online documentation for more information.
setClass("Spectrum",
         representation = representation(
             msLevel="integer",
             peaksCount="integer",
             rt="numeric",
             acquisitionNum="integer",
             scanIndex = "integer",
             tic = "numeric",
             mz = "numeric",
             intensity = "numeric",
             fromFile = "integer",
             centroided = "logical",
             smoothed = "logical",
             polarity="integer",
             "VIRTUAL"),
         contains=c("Versioned"),
         prototype = prototype(
             rt = numeric(),
             polarity = NA_integer_,
             acquisitionNum = NA_integer_,
             msLevel = NA_integer_,
             centroided = NA,
             smoothed = NA,
             peaksCount = 0L,
             tic = 0,
             scanIndex = integer(),
             mz = numeric(),
             intensity = numeric()),
         validity = function(object)
             validSpectrum(object)
         )

setClass("Spectrum2",
         representation = representation(
             merged="numeric",
             precScanNum="integer",
             precursorMz="numeric",
             precursorIntensity = "numeric",
             precursorCharge = "integer",
             collisionEnergy = "numeric"),
         contains=c("Spectrum"),
         prototype = prototype(
             merged = 1,
             acquisitionNum = integer(),
             precScanNum = integer(),
             precursorMz = numeric(),
             precursorIntensity = numeric(),
             msLevel = as.integer(2),
             precursorCharge = integer(),
             collisionEnergy = numeric()),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             msl <- object@msLevel
             if (msl < as.integer(2))
                 msg <- validMsg(msg,
                                 paste0("Object of class ",
                                        class(object),
                                        " but msLevel is ", msl,
                                        " (should be > 1)"))
             if (is.null(msg)) TRUE
             else msg
         })


setClass("Spectrum1",
         contains=c("Spectrum"),
         prototype = prototype(
             polarity=integer(),
             msLevel = as.integer(1)),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             msl <- object@msLevel
             if (msl!=as.integer(1))
                 msg <- validMsg(msg,paste("Object of class",class(object),
                                           "but msLevel is",msl,sep=" "))
             if (is.null(msg)) TRUE
             else msg
         })

##################################################################
## Data Structure for Reporter Ions for labelled MS Quantification
## See online documentation for more information.
setClass("ReporterIons",
         representation = representation(
             name = "character",
             reporterNames = "character",
             description = "character",
             mz = "numeric",
             col = "character",
             width = "numeric"),
         contains = c("Versioned"),
         prototype  =  prototype(
             new("Versioned", versions = c(ReporterIons = "0.1.0")),
             name = character(),
             reporterNames = character(),
             description = character(),
             mz = numeric(),
             col = character(),
             width = numeric()),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (length(object@mz) == 0) {
                 msg <- validMsg(msg, "No reporter ions defined.")
             } else {
                 if (length(object@col) != length(object@mz))
                     warning("Missing color(s) for the reporter ions.")
             }
             if (is.null(msg)) TRUE
             else msg
         })


#####################################################################
## The "MSnSet" Class for MS Proteomics Expression Data and Meta-Data
## See online documentation for more information.
setClass("MSnSet",
         representation = representation(
             experimentData="MIAPE",
             processingData = "MSnProcess",
             qual = "data.frame"),
         contains = "eSet",
         prototype = prototype(
             new("VersionedBiobase",
                 versions=c(classVersion("eSet"),
                     classVersion("pSet"),
                     MSnSet="0.4.0")),
             experimentData=new("MIAPE"),
             annotation="No feature annotation."))

.MSnSetList <-
    setClass("MSnSetList",
             slots = c(x = "list",
                       log = "list",
                       featureData = "DataFrame"),
             contains = "Versioned",
             prototype = prototype(
                 new("Versioned",
                     versions = c(MSnSetList = "0.2.0"))),
             validity = function(object) {
                 msg <- validMsg(NULL, NULL)
                 if (!listOf(object@x, "MSnSet", valid = FALSE))
                     msg <- validMsg(msg, "Not all items are MSnSets.")
                 nvals <- sapply(object@x, validObject)
                 if (!all(nvals))
                     msg <- validMsg(msg,
                                     paste(sum(!nvals),
                                           "MSnSets are not valid."))
                 if (length(object@x) != nrow(object@featureData))
                     msg <- validMsg(msg,
                                     "Data and meta-data dimensions don't match.")
                 if (length(object@x) &&
                     !identical(names(object@x), rownames(object@featureData)))
                     msg <- validMsg(msg,
                                     "Data and meta-data names don't match.")
                 if (is.null(msg)) TRUE
                 else msg
             })

#####################################################################
## Features of interest infrastructure
.FeaturesOfInterest <-
    setClass("FeaturesOfInterest",
             slots = c(
                 description = "character",
                 fnames = "character",
                 date = "character",
                 objpar = "list"),
             contains = "Versioned",
             prototype = prototype(
                 new("Versioned",
                     versions = c(FeaturesOfInterest = "0.1.0"))))

.FoICollection <-
    setClass("FoICollection",
             slots = c(foic = "list"),
             contains = "Versioned",
             prototype = prototype(
                 new("Versioned",
                     versions = c(FeaturesOfInterest = "0.1.0"))))

#####################################################################
## Support for the PSI mzTab format
.MzTab <- setClass("MzTab",
                   slots = c(
                       Metadata = "list",
                       Filename = "character",
                       Proteins = "data.frame",
                       Peptides = "data.frame",
                       PSMs = "data.frame",
                       SmallMolecules = "data.frame",
                       Comments = "character"))

##################################################################
## Container for MSn Experiments Data and Meta-Data enabling/allowing
## to process raw MS files on-the-fly (without the need to keep all
## data in memory).
setClass("OnDiskMSnExp",
         representation=representation(
             spectraProcessingQueue="list",  ## List collecting ProcessingSteps for lazy processing.
             backend="character"  ## That's to eventually add a SQLite backend later...
         ),
         contains=c("MSnExp"),
         prototype = prototype(
             new("VersionedBiobase",
                 versions = c(classVersion("MSnExp"), OnDiskMSnExp="0.0.1")),
             spectraProcessingQueue=list(),
             backend=character()),
         validity=function(object){
             ## Return true if the object is empty.
             if (length(object) == 0)
                 return(TRUE)
             ## Ensure that the files (returned by fileNames) are available
             ## and check also that the featureData contains all the required
             ## information.
             msg <- validMsg(NULL, NULL)
             ## Elements in spectraProcessingQueue have to be ProcessingStep objects.
             if (length(object@spectraProcessingQueue) > 0){
                 isOK <- unlist(lapply(object@spectraProcessingQueue,
                                       function(z)
                                           return(is(z, "ProcessingStep"))))
                 if (any(!isOK))
                     msg <- validMsg(msg,
                                     paste0("Only objects of type 'ProcessingStep'",
                                            " allowed in slot 'spectraProcessingQueue'"))
             }
             ## Check that required columns are present in the featureData:
             msg <- validMsg(msg, validateFeatureDataForOnDiskMSnExp(featureData(object)))
             ## Check if the files do exist.
             theFiles <- fileNames(object)
             for (theF in theFiles){
                 if (!file.exists(theF))
                     msg <- validMsg(msg,
                                     paste0("Required data file '", basename(theF),
                                            "' not found!"))
             }
             ## Some last checks I had to take from the pSet as validObject on OnDiskMSnExp was first
             ## calling the pSet validate method on the MSnExp and that caused a problem since we don't
             ## have an assayData with spectra here (fromFile was trying to get that form there).
             aFileIds <- fromFile(object)
             fFileIds <- fData(object)$fileIdx
             if (length(fFileIds) && any(aFileIds != fFileIds))
                 msg <- validMsg(msg, "Mismatch of files in assayData and processingData.")
             ## Check if the fromFile values match to @files in processingData
             filesProcData <- 1:length(processingData(object)@files)
             if ( !all(unique(sort(aFileIds)) == unique(sort(filesProcData))) )
                 msg <- validMsg(msg, "Spectra file indices in assayData does not match files in processinData.")
             nfilesprocData   <- length(processingData(object)@files)
             nfilesSpectra <- length(unique(aFileIds))
             if (nfilesprocData < nfilesSpectra)
                 msg <- validMsg(msg, "More spectra files in assayData than in processingData.")
             if (length(sampleNames(object)) != nrow(pData(object)))
                 msg <- validMsg(msg, "Different number of samples accoring to sampleNames and pData.")
             ## Check also experimentData:
             if (length(
                 unique(c(length(fileNames(object)),
                          length(experimentData(object)@instrumentManufacturer),
                          length(experimentData(object)@instrumentModel),
                          length(experimentData(object)@ionSource),
                          length(experimentData(object)@analyser),
                          length(experimentData(object)@detectorType)))) != 1)
                 msg <- validMsg(msg, "The number of files does not match the information in experimentData.")
             if (!isOnDisk(object))
                 msg <- validMsg(msg, "Object is not 'onDisk'.")
             if (!isEmpty(object@assayData))
                 msg <- validMsg(msg, "Assaydata is not empty.")
             if (is.null(msg)) {
                 return(TRUE)
             } else {
                 return(msg)
             }
         })

############################################################
## Simple class defining a "processing" step, consisting of a
## function name (@FUN) and optional arguments (@ARG)
setClass("ProcessingStep",
         representation=representation(
             FUN="character",
             ARGS="list"
         ),
         contains="Versioned",
         prototype=prototype(
             new("Versioned",
                 versions=c(ProcessingStep="0.0.1")),
             ARGS=list(),
             FUN=character()
         ),
         validity=function(object){
             msg <- validMsg(NULL, NULL)
             ## Check if function/method exists?
             if(length(object@FUN) > 0){
                 Res <- try(get(object@FUN), silent=TRUE)
                 if(is(Res, "try-error")){
                     msg <- validMsg(msg,
                                     paste0("Function '", object@FUN, "' not found!"))
                 }else{
                     if(!is(Res, "function"))
                         msg <- validMsg(msg,
                                         paste0("'", object@FUN, "' is not a function!"))
                 }
             }
             if(is.null(msg)){
                 return(TRUE)
             }else{
                 return(msg)
             }
         })

#' @title Representation of chromatographic MS data
#'
#' @description The \code{Chromatogram} class is designed to store
#'     chromatographic MS data, i.e. pairs of retention time and intensity
#'     values. Instances of the class can be created with the
#'     \code{Chromatogram} constructor function but in most cases the dedicated
#'     methods for \code{\linkS4class{OnDiskMSnExp}} and
#'     \code{\linkS4class{MSnExp}} objects extracting chromatograms should be
#'     used instead (i.e. the \code{\link{chromatogram}} method).
#'
#' @details The \code{mz}, \code{filterMz}, \code{precursorMz} and
#'     \code{productMz} are stored as a \code{numeric(2)} representing a range
#'     even if the chromatogram was generated for only a single ion (i.e. a
#'     single mz value). Using ranges for \code{mz} values allow this class to
#'     be used also for e.g. total ion chromatograms or base peak chromatograms.
#'
#'     The slots \code{precursorMz} and \code{productMz} allow to represent SRM
#'     (single reaction monitoring) and MRM (multiple SRM) chromatograms. As
#'     example, a \code{Chromatogram} for a SRM transition 273 -> 153 will have
#'     a \code{@precursorMz = c(273, 273)} and a
#'     \code{@productMz = c(153, 153)}.
#'
#' @rdname Chromatogram-class
#'
#' @export
#'
#' @seealso \code{\link{Chromatograms}} for combining \code{Chromatogram} in
#'     a two-dimensional matrix (rows being mz-rt ranges, columns samples).
#'     \code{\link{chromatogram}} for the method to extract chromatogram data
#'     from a \code{\linkS4class{MSnExp}} or \code{\linkS4class{OnDiskMSnExp}}
#'     object.
#'     \code{\link{clean}} for the method to \emph{clean} a \code{Chromatogram}
#'     object.
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Create a simple Chromatogram object.
#' ints <- abs(rnorm(100, sd = 100))
#' rts <- seq_len(length(ints))
#' chr <- Chromatogram(rtime = rts, intensity = ints)
#' chr
#'
#' ## Extract intensities
#' intensity(chr)
#'
#' ## Extract retention times
#' rtime(chr)
#'
#' ## Extract the mz range - is NA for the present example
#' mz(chr)
#'
#' ## plot the Chromatogram
#' plot(chr)
setClass("Chromatogram",
         slots = c(
             rtime = "numeric",
             intensity = "numeric",
             mz = "numeric",
             filterMz = "numeric",
             precursorMz = "numeric", ## Or call that Q1mz?
             productMz = "numeric",   ## Or call that Q3mz?
             fromFile = "integer",
             aggregationFun = "character",
             msLevel = "integer"
         ),
         contains = "Versioned",
         prototype = prototype(
             rtime = numeric(),
             intensity = numeric(),
             mz = c(NA_real_, NA_real_),
             filterMz = c(NA_real_, NA_real_),
             precursorMz = c(NA_real_, NA_real_),
             productMz = c(NA_real_, NA_real_),
             fromFile = integer(),
             aggregationFun = character(),
             msLevel = 1L
         ),
         validity = function(object)
             .validChromatogram(object)
         )

#' @title Container for multiple Chromatogram objects
#'
#' @aliases coerce,matrix,Chromatograms-method
#'
#' @description The \code{Chromatograms} class allows to store
#'     \code{\link{Chromatogram}} objects in a \code{matrix}-like
#'     two-dimensional structure.
#'
#' @details The \code{Chromatograms} class extends the base \code{matrix} class
#'     and hence allows to store \code{\link{Chromatogram}} objects in a
#'     two-dimensional array. Each row is supposed to contain
#'     \code{Chromatogram} objects for one MS data \emph{slice} with a common
#'     m/z and rt range. Columns contain \code{Chromatogram} objects from the
#'     same sample.
#'
#' @export
#'
#' @rdname Chromatograms-class
#'
#' @seealso \code{\link{Chromatogram}} for the class representing chromatogram
#'     data.
#'     \code{\link{chromatogram}} for the method to extract a
#'     \code{Chromatograms} object from a \code{\linkS4class{MSnExp}} or
#'     \code{\linkS4class{OnDiskMSnExp}} object.
#'     \code{\link{readSRMData}} for the function to read chromatographic data
#'     of an SRM/MRM experiment.
#'
#' @author Johannes Rainer
#'
#' @examples
#' ## Creating some chromatogram objects to put them into a Chromatograms object
#' ints <- abs(rnorm(25, sd = 200))
#' ch1 <- Chromatogram(rtime = 1:length(ints), ints)
#' ints <- abs(rnorm(32, sd = 90))
#' ch2 <- Chromatogram(rtime = 1:length(ints), ints)
#' ints <- abs(rnorm(19, sd = 120))
#' ch3 <- Chromatogram(rtime = 1:length(ints), ints)
#' ints <- abs(rnorm(21, sd = 40))
#' ch4 <- Chromatogram(rtime = 1:length(ints), ints)
#'
#' ## Create a Chromatograms object with 2 rows and 2 columns
#' chrs <- Chromatograms(list(ch1, ch2, ch3, ch4), nrow = 2)
#' chrs
#'
#' ## Extract the first element from the second column. Extracting a single
#' ## element always returns a Chromatogram object.
#' chrs[1, 2]
#'
#' ## Extract the second row. Extracting a row or column (i.e. multiple elements
#' ## returns by default a list of Chromatogram objects.
#' chrs[2, ]
#'
#' ## Extract the second row with drop = FALSE, i.e. return a Chromatograms
#' ## object.
#' chrs[2, , drop = FALSE]
#'
#' ## Replace the first element.
#' chrs[1, 1] <- ch3
#' chrs
#'
#' ## Add a pheno data.
#' pd <- data.frame(name = c("first sample", "second sample"),
#'     idx = 1:2)
#' pData(chrs) <- pd
#'
#' ## Column names correspond to the row names of the pheno data
#' chrs
#'
#' ## Access a column within the pheno data
#' chrs$name
#'
#' ## Access the m/z ratio for each row; this will be NA for the present
#' ## object
#' mz(chrs)
setClass("Chromatograms",
         contains = "matrix",
         slots = c(phenoData = "NAnnotatedDataFrame",
                   featureData = "AnnotatedDataFrame"),
         prototype = prototype(
             matrix(ncol = 0, nrow = 0),
             phenoData = new("NAnnotatedDataFrame",
                             dimLabels = c("sampleNames", "sampleColumns")),
             featureData = new("AnnotatedDataFrame",
                               dimLabels = c("featureNames", "featureColumns"))
         ),
         validity = function(object)
             .validChromatograms(object)
         )

#' @name Spectra
#'
#' @aliases Spectra-class show,Spectra-method coerce,Spectra,list-method coerce,Spectra,MSnExp-method
#'
#' @title List of Spectrum objects along with annotations
#'
#' @description
#'
#' `Spectra` objects allow to collect one or more [Spectrum-class] object(s)
#' ([Spectrum1-class] or [Spectrum2-class]) in a `list`-like structure with
#' the possibility to add arbitrary annotations to each individual
#' `Spectrum` object. These can be accessed/set with the [mcols()] method.
#'
#' `Spectra` objects can be created with the `Spectra` function.
#'
#' Functions to access the individual spectra's attributes are available
#' (listed below).
#'
#' @details
#'
#' `Spectra` inherits all methods from the [SimpleList] class of the
#' `S4Vectors` package. This includes `lapply` and other data manipulation
#' and subsetting operations.
#'
#' @param object For all functions: a `Spectra` object.
#'
#' @param x For all functions: a `Spectra` object.
#'
#' @md
#'
#' @rdname Spectra
NULL

.Spectra <- setClass("Spectra",
                     contains = "SimpleList",
                     prototype = prototype(elementType = "Spectrum")
                     )

setValidity("Spectra", function(object) {
    ## All elements in the list have to be Spectrum objects.
    msg <- character()
    if (any(vapply(object, function(z) !is(z, "Spectrum"), logical(1))))
        msg <- c(msg, "All elements have to be Spectrum objects")
    if (length(msg)) msg else TRUE
})
