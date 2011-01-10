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
           centroided="logical",
           normalised="logical",
           xcmsVersion="character",
           MSnbaseVersion="character"),
         contains=c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(MSnProcess="0.1.0")),
           processing=character(),
           removedPeaks=character(),
           xcmsVersion=as.character(packageDescription("xcms",fields="Version")),
           ## will have to check whether this is a problem during 
           ## package building, checking when packahe not yet installed
           ## as well as during first installation
           ## Put this in an initialise() method!
           MSnbaseVersion=ifelse(
             is.na(packageDescription("MSnbase",fields="Version")),
             "0.0.0",
             as.character(packageDescription("MSnbase",fields="Version")))
           )
         )

#################################################################
## The 'Minimum Information About a Proteomics Experiment' Class
## See online documentation for more information.
setClass("MIAPE",
         representation=representation(
           proteomicsData="character"),
         contains=c("MIAME"),
         prototype = prototype(
           new("Versioned", versions=c(MIAPE="0.0.1")),
           proteomicsData="Will contain more MIAPE data."
           )
         )

############################################################################
## NAnnotatedDataFrame: As Biobase's AnnotatedDataFrame, it is composed of
## a data.frame, with annotations about columns named
## in the data slot contained in the metadata slot.
## In addition, it contains a multiplex slot to make explicite that
## the AnnotatedDataFrame is applied to a set of mulitplexed tags.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setClass("NAnnotatedDataFrame",
         representation(multiplex="character"),
         contains = c("AnnotatedDataFrame"),
         prototype = prototype(
           new("Versioned", versions=list(NAnnotatedDataFrame="0.0.1"))))


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
                        protocolData = "AnnotatedDataFrame",
                        process = "MSnProcess",
                        "VIRTUAL"),
         contains = "VersionedBiobase",
         prototype = prototype(
           new("VersionedBiobase", versions=c(pSet="0.0.1")),
           assayData = new.env(),
           phenoData = new("NAnnotatedDataFrame",
             dimLabels=c("sampleNames", "sampleColumns")), ## here, it's rather files than samples...
           featureData = new("AnnotatedDataFrame",
             dimLabels=c("featureNames", "featureColumns")),
           protocolData = new("AnnotatedDataFrame",
             dimLabels=c("sampleNames", "sampleColumns"))))

##################################################################
## Container for MSn Experiments Data and Meta-Data
## See online documentation for more information.
setClass("MSnExp",
         contains=c("pSet"),
         prototype = prototype(
           new("VersionedBiobase",
               versions=c(classVersion("pSet"), MSnExp="0.3.0")),
           experimentData=new("MIAPE")))


## setClass("MSnExp",
##          representation = representation(
##            spectra="list",
##            process="MSnProcess",
##            fromFile="numeric",
##            files="character"),
##          contains=c("eSet"),
##          prototype = prototype(
##            spectra=list(),
##            process=new("MSnProcess"),
##            experimentData=new("MIAPE"),
##            fromFile=numeric(),
##            files=character(),           
##            new("VersionedBiobase",
##                versions=c(classVersion("eSet"), MSnExp="0.2.0")))
##          )


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
           fromFile = "integer", ## added to v0.1.1 to replace fromFile in MSnExp
           "VIRTUAL"),
         contains=c("Versioned"),
         prototype = prototype(
           new("Versioned", versions=c(Spectrum="0.1.1")),
           rt = numeric(),
           acquisitionNum = integer(),
           msLevel = integer(),
           peaksCount = integer(), 
           scanindex = integer(),
           mz = numeric(),
           intensity = numeric()),
         validity = function(object) {
           msg <- validMsg(NULL, NULL)
           if (any(is.na(object@intensity)))
             msg <- validMsg(msg,"'NA' MZ values found.")
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
           process = "MSnProcess",
           qual = "data.frame"),
         contains = c("ExpressionSet"),
         prototype = prototype(
           new("VersionedBiobase",
               versions=c(classVersion("ExpressionSet"),classVersion("MSnSet"), MSnSet="0.3.0")),
           experimentData=new("MIAPE")))

