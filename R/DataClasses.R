#################################################################
## The 'Minimum Information About a Proteomics Experiment' Class
## See online documentation for more information.
setClass("MIAPE",
         representation=representation(description="character"),
         contains=c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(MIAPE="0.0.1")),
           description="Will contain MIAPE data."
           )
         )

######################################################################
## MSnProcess: Container for MSnExp and MSnSet processing information
## See online documentation for more information.
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
           MSnbaseVersion=ifelse(
             all(is.na(packageDescription("MSnbase"))),
             "0.0.0",
             as.character(packageVersion("MSnbase")))
           )
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
           new("Versioned",
               versions=c(classVersion("Spectrum"), Spectrum2="0.1.0")),
           merged = 1,
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
## Container for MSn Experiments Data and Meta-Data
## See online documentation for more information.
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
           spectra=list(),
           process=new("MSnProcess"),
           proteomicsData=new("MIAPE"),
           description="Provide a short description of your experiment here.",
           fromFile=numeric(),
           files=character(),           
           new("VersionedBiobase",
               versions=c(classVersion("eSet"), MSnExp="0.2.0")))
         )

 
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
               warning("Missing colors for the reporter ions.")
           }
           if (is.null(msg)) TRUE
           else msg
         })



#####################################################################
## The "MSnSet" Class for MS Proteomics Expression Data and Meta-Data
## See online documentation for more information.
setClass("MSnSet",
         representation = representation(
           qual="data.frame",
           process="MSnProcess",
           proteomicsData="MIAPE",
           description="character",
           files="character"),
         contains = c("ExpressionSet"),
         prototype = prototype(
           proteomicsData=new("MIAPE"),
           process=new("MSnProcess"),
           new("VersionedBiobase",
               versions=c(classVersion("ExpressionSet"), MSnSet="0.2.0")))
         )


