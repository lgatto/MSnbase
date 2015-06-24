######################################################################
## MSnProcess: Container for MSnExp and MSnSet processing information
## See online documentation for more information.
setClass("MSnProcess",
         representation = representation(
             files="character",
             processing="character",
             merged="logical",
             cleaned="logical",
             removedPeaks="character",
             smoothed="logical",
             trimmed="numeric",
             normalised="logical",
             MSnbaseVersion="character"),
         contains=c("Versioned"),
         prototype = prototype(
             new("Versioned", versions=c(MSnProcess="0.1.3")),
             processing=character(),
             files=character(),
             trimmed=numeric(),
             removedPeaks=character(),
             MSnbaseVersion=character()) ## set in initialize()
         )

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
##########################
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
             ## 3. Post-source componentry           
             analyser = "character", ## Quad, TOF, Trap, ... 
             analyserDetails = "character",
             ## 3. Post-source componentry - (d) Collision cell
             collisionGas = "character",
             collisionPressure = "numeric",
             collisionEnergy = "character",
             ## 3. Post-source component — (f) Detectors
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
             new("Versioned", versions=c(classVersion("MIAxE"), MIAPE="0.2.2")),
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
             other = list())
         )

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
             new("Versioned", versions=list(NAnnotatedDataFrame="0.0.3")),
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
                        .cache = "environment", ## locked
                        "VIRTUAL"),
         contains = "Versioned",
         prototype = prototype(
             new("Versioned", versions=c(pSet="0.1.0")),
             assayData = new.env(parent=emptyenv()), 
             experimentData = new("MIAPE"),
             phenoData = new("NAnnotatedDataFrame",
                 dimLabels=c("sampleNames", "fileNumbers")), 
             featureData = new("AnnotatedDataFrame",
                 dimLabels=c("featureNames", "featureColumns")),
             protocolData = new("AnnotatedDataFrame",
                 dimLabels=c("sampleNames", "sampleColumns")))
         )

##################################################################
## Container for MSn Experiments Data and Meta-Data
## See online documentation for more information.
setClass("MSnExp",
         contains=c("pSet"),
         prototype = prototype(
             new("VersionedBiobase",
                 versions=c(classVersion("pSet"), MSnExp="0.3.0")),
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
             "VIRTUAL"),
         contains=c("Versioned"),
         prototype = prototype(
             new("Versioned", versions=c(Spectrum="0.2.0")),
             rt = numeric(),
             acquisitionNum = integer(),
             msLevel = integer(),
             centroided = FALSE,
             peaksCount = 0L,
             tic = 0L,
             scanIndex = integer(),
             mz = numeric(),
             intensity = numeric()),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (any(is.na(intensity(object))))
                 msg <- validMsg(msg,"'NA' intensities found.")
             if (any(is.na(mz(object))))
                 msg <- validMsg(msg,"'NA' M/Z found.")
             if (any(intensity(object) < 0))
                 msg <- validMsg(msg,"Negative intensities found.")
             if (any(mz(object)<0))
                 msg <- validMsg(msg,"Negative M/Z found.")
             if (length(object@mz) != length(object@intensity))
                 msg <- validMsg(msg,"Unequal number of MZ and intensity values.")
             if (length(mz(object)) != peaksCount(object))
                 msg <- validMsg(msg,"Peaks count does not match up with number of MZ values.")
             if (is.null(msg)) TRUE
             else msg
         })

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
             new("Versioned",
                 versions=c(classVersion("Spectrum"), Spectrum2="0.1.0")),
             merged = 1,
             centroided = FALSE,
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
             if (msl!=as.integer(2))
                 msg <- validMsg(msg,paste("Object of class",class(object),
                                           "but msLevel is",msl,sep=" "))
             if (is.null(msg)) TRUE
             else msg
         })


setClass("Spectrum1",
         representation = representation(polarity="integer"),
         contains=c("Spectrum"),
         prototype = prototype(
             new("Versioned", versions=c(classVersion("Spectrum"), Spectrum1="0.1.0")),
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
             name="character",
             reporterNames="character",
             description="character",
             mz="numeric",
             col="character",
             width="numeric"),
         contains=c("Versioned"),
         prototype = prototype(
             new("Versioned", versions=c(ReporterIons="0.1.0")),
             name=character(),
             reporterNames=character(),
             description=character(),
             mz=numeric(),
             col=character(),
             width=numeric()),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (length(object@mz)==0) {
                 msg <- validMsg(msg,"No reporter ions defined.")
             } else {
                 if (length(object@col)!=length(object@mz))
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
                       log = "list"),
             contains = "Versioned",
             prototype = prototype(
                 new("Versioned",
                     versions = c(MSnSetList = "0.1.0"))),
             validity = function(object) {
                 msg <- validMsg(NULL, NULL)
                 if (!listOf(object@x, "MSnSet"))
                     msg <- validMsg(msg, "Not all items are MSnSets.")
                 nvals <- sapply(object@x, validObject)
                 if (!all(nvals))
                     msg <- validMsg(msg,
                                     paste(sum(!nvals),
                                           "MSnSets are not valid."))
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
