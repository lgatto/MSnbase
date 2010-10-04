##' The Minimum Information About a Proteomics Experiment
##'
##' The current dummy implementation is composed of a free-text
##' description of the experiment. The class will be updated
##' at a later stage to be MIAPE compliant.
##'
##' @title The MIAPE Class
##' @name MIAPE
##' @slot description a character string describing the experiment
##' @aliases MIAPE-class
##' @aliases MIAPE
##' @references
##' About MIAPE: \url{http://www.psidev.info/index.php?q=node/91}
##' Reporting guideline: \url{http://www.psidev.info/index.php?q=node/60}
##' @author Laurent Gatto <lg390@@cam.ac.uk>
##' @examples
##' new("MIAPE",description="Quantitative MSMS experiment using iTRAQ 4-plex on QTop Premier.")
##' @TODO
##' This is a temporary implementation, that will most likely
##' be updated using the AnnotedDataFrame class. MIAPE is currently 
##' just a checklist of items, but may be more formally described later.
##' @keywords classes
##' @docType class 
##' @exportClass MIAPE
setClass("MIAPE",
         representation=representation(description="character"),
         contains=c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(MIAPE="0.0.1")),
           description="Will contain MIAPE data."
           )
         )


##' Container for MS Spectra
##'
##' Virtual container for spectrum data common
##' to all different types of spectra. A \code{Spectrum}
##' object can not be directly instanciated. Use
##' \code{"\linkS4class{Spectrum1}"} and
##' \code{"\linkS4class{Spectrum2}"} instead. 
##'
##' @title The Virtual Spectrum Class
##' @name Spectrum
##' @slot msLevel an integer indication the MS level; 
##' 1 for MS1 level \code{Spectrum1} objects and
##' 2 for MSMSM \code{Spectrum2} objects.
##' @slot peaksCount an integer indicating the number
##' of MZ peaks.
##' @slot rt a numeric value indicating the retention
##' time (in seconds) for the current ions.
##' @slot acquisitionNum an integer corresponding to
##' the acquisition number of the current spectrum.
##' @slot scanIndex an integer indicating the scan
##' index of the current spectrum.
##' @slot mz a numeric of length equal to the peaks count
##' indicating the MZ values that have been measured for
##' the current ion.
##' @slot intensity a numeric of same length as \code{mz}
##' indicating the intensity at which each \code{mz} datum
##' has been measured.
##' @aliases Spectrum-class
##' @aliases Spectrum
##' @author Laurent Gatto <lg390@@cam.ac.uk>
##' @note This is a virtual class and can not be
##' instanciated directly
##' @seealso Instaciable sub-classes
##' \code{"\linkS4class{Spectrum1}"}
##' and \code{"\linkS4class{Spectrum2}"}
##' for MS1 and MS2 spectra.
##' @keywords classes
##' @docType class 
##' @exportClass Spectrum
setClass("Spectrum",
         representation = representation(
           msLevel="integer",
           peaksCount="integer",
           rt="numeric",
           acquisitionNum="integer",
           scanIndex = "integer",
           mz = "numeric",
           intensity = "numeric",
           "VIRTUAL"),
         contains=c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(Spectrum="0.1.0")),
           rt = numeric(),
           acquisitionNum = integer(),
           msLevel = integer(),
           peaksCount = integer(), 
           scanindex = integer(),
           mz = numeric(),
           intensity = numeric()),
         validity = function(object) {
           msg <- validMsg(NULL, NULL)
           if (length(object@mz)!=length(object@intensity))
             msg <- validMsg(msg,"Unequal number of MZ and intensity values.")
           if (length(object@mz)!=object@peaksCount)
             msg <- validMsg(msg,"Peaks count does not match up with number of MZ values.")           
           msl <- object@msLevel
           if (length(unique(msl))!=1) 
             warning(paste("Different MS levels in ",class(object),
                           " object:",unique(msl)))
           if (is.null(msg)) TRUE
           else msg
         })

setClass("Spectrum2",
         representation = representation(
           merged="numeric",
           ms1scan="integer",
           precursorMz="numeric",
           precursorIntensity = "numeric",
           precursorCharge = "integer",
           scanindex = "integer",
           collisionEnergy = "numeric"),
         contains=c("Spectrum"),
         prototype = prototype(
           new("Versioned", versions=c(Spectrum1="0.1.0")),
           merged = numeric(),
           acquisitionNum = integer(),
           ms1scan = integer(),
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
           new("Versioned", versions=c(Spectrum1="0.1.0")),
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


setClass("MSnProcess",
         representation = representation(
           processing="character",
           merged="logical",
           cleaned="logical",
           removedPeaks="character",
           smoothed="logical",
           centroided="logical",
           normalised="logical",
           xcmsVersion="character",
           MSnbaseVersion="character"),
         contains=c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(MSnProcess="0.1.0")),
           processing=character(),
           merged=FALSE,
           cleaned=FALSE,
           removedPeaks=character(),
           smoothed=FALSE,
           centroided=FALSE,
           normalised=FALSE,
           xcmsVersion=as.character(packageVersion("xcms")),
           ## will have to check whether this is a problem during 
           ## package building, checking when packahe not yet installed
           ## as well as during first installation
           MSnbaseVersion=ifelse(is.na(packageDescription("MSnbase")),"0.0.0",
             as.character(packageVersion("MSnbase")))
           )
         )
         

setClass("MSnExp",
         representation = representation(
           spectra="list",
           process="MSnProcess",
           proteomicsData="MIAPE",
           description="character",
           fromFile="numeric",
           files="character"),
         contains=c("eSet"),
         prototype = prototype(
           new("VersionedBiobase",versions=c(MSnExp="0.2.0")),
           description="Short description of the data.",
           fromFile=numeric(),
           files=character())
         )

setClass("ReporterIons",
         representation = representation(
           name="character",
           description="character",
           mz="numeric",
           col="character",
           width="numeric"),
         contains=c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(ReporterIons="0.1.0")),
           name=character(),
           description=character(),
           mz=numeric(),
           col=character(),
           width=numeric()),
         validity = function(object) {
           msg <- validMsg(NULL, NULL)
           if (length(object@mz)==0)
             msg <- validMsg(msg,"No reporter ions defined.")
           if (length(object@mz)>0) {
             if (length(object@col)!=length(object@mz))
               warning("Missing colors for the reporter ions.")
           }
           if (is.null(msg)) TRUE
           else msg
         })


iTRAQ4 <- new("ReporterIons",
              description="4-plex iTRAQ",
              name="iTRAQ4",
              mz=c(114.13,115.13,116.13,117.13),
              col=c("red","green","blue","yellow"),
              width=0.05)

iTRAQ5 <- new("ReporterIons",
              description="4-plex iTRAQ with isobaric tag",
              name="iTRAQ4",
              mz=c(114.13,115.13,116.13,117.13,145.13),
              col=c("red","green","blue","yellow","grey"),
              width=0.05)

TMT6 <- new("ReporterIons",
            description="6-plex TMT tags",
            name="TMT6",
            mz=c(126.22,127.21,128.20,129.20,130.19,131.18),
            col=c("red","purple","blue","steelblue","green","yellow"),
            width=0.05)

TMT7 <- new("ReporterIons",
            description="6-plex TMT tags with isobaric tag",
            name="TMT7",
            mz=c(126.22,127.21,128.20,129.20,130.19,131.18,229.26),
            col=c("red","purple","blue","steelblue","green","yellow","grey"),
            width=0.05)


setClass("MSnQual",
         representation = representation(
           qc="data.frame"
           ## add metadata and process
           ),
         contains=c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(MSnQual="0.1.0")))
         )


setClass("MSnSet",
         representation = representation(
           proteomicsData="MIAPE",
           process="MSnProcess"),
         contains = c("ExpressionSet"),
         prototype = prototype(
           new("VersionedBiobase",versions=c(MSnSet="0.1.0")))
         )




