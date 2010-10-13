##################################################################
## Methods for MIAPE class
setMethod("show","MIAPE",function(object) show.MIAPE(object))


##################################################################
## Methods for ReporterIons class
setMethod("show","ReporterIons",function(object) show.ReporterIons(object))
setMethod("[","ReporterIons",
          function(x,i,j="missing",drop="missing") "[.ReporterIons"(x,i))
setMethod("length","ReporterIons",function(x) length(x@mz))
setMethod("reporterNames","ReporterIons", function(object) object@reporterNames)
setMethod("reporterNames<-","ReporterIons",
          function(object,value="list") object@reporterNames <- value)
setReplaceMethod("reporterNames",
                 signature(object="ReporterIons",
                           value="character"),
                 function(object, value) {
                   if (length(value)!=length(object))
                     stop(paste("Please provide names for",
                                length(object),
                                "reporters",sep=" "))
                   object@reporterNames = value
                   if (validObject(object))
                     return(object)
                 })
setMethod("initialize", "ReporterIons",
          function(.Object,...) {
            .Object <- callNextMethod()
            if (length(.Object@mz)!=length(.Object@reporterNames)) {
              warning("Setting reporter name.")
              .Object@reporterNames <- paste(.Object@name,.Object@mz,sep=".")
            }
            .Object
          })


##################################################################
## Methods for Spectrum class and children
setMethod("show","Spectrum",function(object) show.Spectrum(object))
setMethod("plot",c("Spectrum","missing"),
          function(x,y,...) {
            if (msLevel(x)==1) plot.Spectrum1(x,...)
            else plot.Spectrum2(x,...)
          })
setMethod("clean","Spectrum",function(object) clean.Spectrum(object))
setMethod("removePeaks","Spectrum",
          function(object,t) removePeaks.Spectrum(object,t))
setMethod("bg.correct","Spectrum",
          function(object,bg) bg.correct.Spectrum(object,bg=-1))
setMethod("precursorMz","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@precursorMz)
            stop("No precursor MZ value for MS1 spectra.")
          })
setMethod("precursorCharge","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@precursorCharge)
            stop("No precursor charge value for MS1 spectra.")
          })
setMethod("acquisitionNum","Spectrum",function(object) object@acquisitionNum)
## setMethod("ms1scanNum","Spectrum",
##           function(object) {
##             if (msLevel(object)>1) 
##               return(object@ms1scan)
##             stop("This is already an MS1 spectrum.")
##           })
setMethod("rtime","Spectrum",function(object) object@rt)
setMethod("peaksCount","Spectrum",function(object) object@peaksCount)
setMethod("msLevel","Spectrum",function(object) object@msLevel)
setMethod("collisionEnergy","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@collisionEnergy)
            stop("No collision energy for MS1 spectra.")
          })
setMethod("intensity","Spectrum",function(object) object@intensity)
setMethod("mz","Spectrum",function(object) object@mz)
setMethod("trimMz","Spectrum",
          function(object,mzlim,...) trimMz.Spectrum(object,mzlim))
setMethod("quantify","Spectrum",
          function(object,reporters,
                   method=c("trapezoidation","max","sum")) 
          quantify.Spectrum(object,reporters,match.arg(method)))
setMethod("curveStats","Spectrum",
          function(object,reporters) curveStats.Spectrum(object,reporters))
setMethod("precursorCharge<-","Spectrum2",
          function(object,value="integer") object@precursorCharge <- value)
setReplaceMethod("precursorCharge",
                 signature(object="Spectrum2",
                           value="integer"),
                 function(object, value) {
                   object@precursorCharge = value
                   if (validObject(object))
                     return(object)
                 })

##################################################################
## Methods for MSnExp class
setMethod("show","MSnExp",function(object) show.MSnExp(object))
setMethod("header","MSnExp",function(object) header.MSnExp(object))
setMethod("length","MSnExp",function(x) length(spectra(x)))
setMethod("plot",c("MSnExp","missing"),
          function(x,y,...) plot.MSnExp(x,...))
setMethod("clean","MSnExp",function(object) clean.MSnExp(object))
setMethod("spectra","MSnExp",function(object) object@spectra)
setMethod("spectra<-","MSnExp",
          function(object,value="list") object@spectra <- value)
setReplaceMethod("spectra",
                 signature(object="MSnExp",
                           value="list"),
                 function(object, value) {
                   object@spectra = value
                   if (validObject(object))
                     return(object)
                 })
setMethod("[","MSnExp",
          function(x,i,j="missing",drop="missing") "[.MSnExp"(x,i))
setMethod("removePeaks","MSnExp",
          function(object,t,verbose) removePeaks.MSnExp(object,t,verbose=TRUE))
setMethod("bg.correct","MSnExp",
          function(object,bg,verbose)
          bg.correct.MSnExp(object,bg=-1,verbose=TRUE))
setMethod("precursorMz","MSnExp",
          function(object) {
            ## this assumes that if first spectrum
            ## has msLevel>1, all have
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), precursorMz))
            stop("No precursor MZ value for MS1 spectra.")
          })
setMethod("precursorCharge","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), precursorCharge))
            stop("No precursor MZ value for MS1 spectra.")
          })
setMethod("acquisitionNum","MSnExp",
          function(object) sapply(spectra(object), acquisitionNum))
## setMethod("ms1scanIdx","MSnExp",
##           function(object) {
##             if (msLevel(object)[1]>1) 
##               return(sapply(spectra(object), ms1scanNum))
##             stop("This experiment contains MS1 spectra.")
##           })
setMethod("rtime","MSnExp",function(object) sapply(spectra(object),rtime))
setMethod("peaksCount","MSnExp",
          function(object) sapply(spectra(object),peaksCount))
setMethod("msLevel","MSnExp",
          function(object) sapply(spectra(object),msLevel))
setMethod("collisionEnergy","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object),collisionEnergy))
            stop("No collision energy for MS1 spectra.")
          })
setMethod("intensity","MSnExp",
          function(object) lapply(object@spectra,intensity))
setMethod("mz","MSnExp",function(object) lapply(object@spectra,mz))
setMethod("trimMz","MSnExp",
          function(object,mzlim,...) {
            object@spectra <- lapply(spectra(object),trimMz,mzlim,...)
            return(object)
          })
setMethod("quantify","MSnExp",
          function(object,reporters,
                   method=c("trapezoidation","max","sum"),
                   verbose=TRUE) {
            ## this assumes that if first spectrum
            ## has msLevel>1, all have
            if (msLevel(object)[1]<2) 
              stop("No quantification for MS1 data implemented.")
            quantify.MSnExp(object,reporters,match.arg(method),verbose)})
setMethod("curveStats","MSnExp",
          function(object,reporters,verbose=TRUE) {
            ifelse(verbose,progress <- "text",progress <- "none")
            l <- llply(object@spectra,curveStats,reporters,.progress=progress)
            qdfr <- l[[1]]
            for (i in 2:length(l)) 
              qdfr <- rbind(qdfr,l[[i]])
            return(qdfr)
          })
##################################################################
## Methods for MSnProcess class
setMethod("show","MSnProcess",function(object) show.MSnProcess(object))


##################################################################
## Methods for MSnSet class
setMethod("normalise","MSnSet",
          function(object,method=c("sum","max"))
          normalise.MSnSet(object,match.arg(method))
          )
setMethod("dim","MSnSet",function(x) dim(exprs(x)))
setMethod("ratios","MSnSet",function(object) ratios.MSnSet(object))
setMethod("qual","MSnSet", function(object) object@qual)
## Not sure about these...
setMethod("featureNames<-","MSnSet",
          function(object,value="character") object@features <- value)
setReplaceMethod("featureNames",
                 signature(object="MSnSet",
                           value="character"),
                 function(object, value) {
                   object@features = value
                   if (validObject(object))
                     return(object)
                 })


