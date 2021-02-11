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
                        phenoData = "AnnotatedDataFrame",
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
             phenoData = new("AnnotatedDataFrame",
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
                       MoleculeFeatures = "data.frame",
                       MoleculeEvidence = "data.frame",
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

setClass("MChromatograms",
         contains = "matrix",
         slots = c(phenoData = "AnnotatedDataFrame",
                   featureData = "AnnotatedDataFrame"),
         prototype = prototype(
             matrix(ncol = 0, nrow = 0),
             phenoData = new("AnnotatedDataFrame",
                             dimLabels = c("sampleNames", "sampleColumns")),
             featureData = new("AnnotatedDataFrame",
                               dimLabels = c("featureNames", "featureColumns"))
         ),
         validity = function(object)
             .validMChromatograms(object)
         )

#' @name MSpectra
#'
#' @aliases MSpectra-class show,MSpectra-method coerce,MSpectra,list-method coerce,MSpectra,MSnExp-method
#'
#' @title List of Spectrum objects along with annotations
#'
#' @description
#'
#' `MSpectra` (Mass Spectra) objects allow to collect one or more
#' [Spectrum-class] object(s) ([Spectrum1-class] or [Spectrum2-class]) in
#' a `list`-like structure with the possibility to add arbitrary annotations
#' to each individual `Spectrum` object. These can be accessed/set with
#' the [mcols()] method.
#'
#' `MSpectra` objects can be created with the `MSpectra` function.
#'
#' Functions to access the individual spectra's attributes are available
#' (listed below).
#'
#' @details
#'
#' `MSpectra` inherits all methods from the [SimpleList] class of the
#' `S4Vectors` package. This includes `lapply` and other data manipulation
#' and subsetting operations.
#'
#' @param object For all functions: a `MSpectra` object.
#'
#' @param x For all functions: a `MSpectra` object.
#'
#' @md
#'
#' @rdname MSpectra
NULL

.MSpectra <- setClass("MSpectra",
                     contains = "SimpleList",
                     prototype = prototype(elementType = "Spectrum")
                     )

setValidity("MSpectra", function(object) {
    ## All elements in the list have to be Spectrum objects.
    msg <- character()
    if (any(vapply(object, function(z) !is(z, "Spectrum"), logical(1))))
        msg <- c(msg, "All elements have to be Spectrum objects")
    if (length(msg)) msg else TRUE
})
