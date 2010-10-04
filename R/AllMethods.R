##' @nord
setMethod("show","MIAPE",function(object) show.MIAPE(object))
##' @nord
setMethod("show","ReporterIons",function(object) show.ReporterIons(object))
##' @nord
setMethod("show","Spectrum",function(object) show.Spectrum(object))
##' @nord
setMethod("show","MSnExp",function(object) show.MSnExp(object))
##' @nord
setMethod("show","MSnProcess",function(object) show.MSnProcess(object))
##' @nord
setMethod("show","MSnQual",function(object) show.MSnQual(object))
##' @nord
setMethod("header","MSnExp",function(object) header.MSnExp(object))
##' @nord
setMethod("length","MSnExp",function(x) length(spectra(x)))

##' @nord
setMethod("plot",c("Spectrum","missing"),
          function(x,y,...) {
            if (msLevel(x)==1) plot.Spectrum1(x,...)
            else plot.Spectrum2(x,...)
          })

##' @nord
setMethod("plot",c("MSnExp","missing"),
          function(x,y,...) plot.MSnExp(x,...))


##' @nord
setMethod("clean","MSnExp",function(object) clean.MSnExp(object))
##' @nord
setMethod("clean","Spectrum",function(object) clean.Spectrum(object))
##' @nord
setMethod("spectra","MSnExp",function(object) object@spectra)
##' @nord
setMethod("spectra<-","MSnExp",
          function(object,value="list") object@spectra <- value)
##' @nord
setReplaceMethod("spectra",
                 signature(object="MSnExp",
                           value="list"),
                 function(object, value) {
                   object@spectra = value
                   if (validObject(object))
                     return(object)
                 })
##' @nord
setMethod("[","MSnExp",function(x,i,j="missing",drop="missing") "[.MSnExp"(x,i))

##' @nord
setMethod("removePeaks","Spectrum",
          function(object,t) removePeaks.Spectrum(object,t))
##' @nord
setMethod("removePeaks","MSnExp",
          function(object,t,verbose) removePeaks.MSnExp(object,t,verbose=TRUE))
##' @nord
setMethod("bg.correct","Spectrum",
          function(object,bg) bg.correct.Spectrum(object,bg="min"))
##' @nord
setMethod("bg.correct","MSnExp",
          function(object,bg,verbose) bg.correct.MSnExp(object,bg="min",verbose=TRUE))


##' @nord
setMethod("precursorMz","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@precursorMz)
            stop("No precursor MZ value for MS1 spectra.")
          })

##' @nord
setMethod("precursorMz","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), precursorMz))
            stop("No precursor MZ value for MS1 spectra.")
          })

##' @nord
setMethod("precursorCharge","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@precursorCharge)
            stop("No precursor charge value for MS1 spectra.")
          })


##' @nord
setMethod("precursorCharge","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), precursorCharge))
            stop("No precursor MZ value for MS1 spectra.")
          })

##' @nord
setMethod("acquisitionNum","Spectrum",function(object) object@acquisitionNum)
##' @nord
setMethod("acquisitionNum","MSnExp",
          function(object) sapply(spectra(object), acquisitionNum))

##' @nord
setMethod("ms1scanNum","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@ms1scan)
            stop("This is already an MS1 spectrum.")
          })

##' @nord
setMethod("ms1scanNum","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), ms1scanNum))
            stop("This experiment contains MS1 spectra.")
          })

##' @nord
setMethod("rtime","Spectrum",function(object) object@rt)
##' @nord
setMethod("rtime","MSnExp",function(object) sapply(spectra(object),rtime))

##' @nord
setMethod("peaksCount","Spectrum",function(object) object@peaksCount)
##' @nord
setMethod("peaksCount","MSnExp",function(object) sapply(spectra(object),peaksCount))

##' @nord
setMethod("msLevel","Spectrum",function(object) object@msLevel)
##' @nord
setMethod("msLevel","MSnExp",function(object) sapply(spectra(object),msLevel))

##' @nord
setMethod("collisionEnergy","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@collisionEnergy)
            stop("No collision energy for MS1 spectra.")
          })
          
##' @nord
setMethod("collisionEnergy","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object),collisionEnergy))
            stop("No collision energy for MS1 spectra.")
          })

##' @nord
setMethod("[","ReporterIons",function(x,i,j="missing",drop="missing") "[.ReporterIons"(x,i))
##' @nord
setMethod("length","ReporterIons",function(x) length(x@mz))

##' @nord
setMethod("intensity","Spectrum",function(object) object@intensity)
##' @nord
setMethod("intensity","MSnExp",function(object) lapply(object@spectra,intensity))
##' @nord
setMethod("mz","Spectrum",function(object) object@mz)
##' @nord
setMethod("mz","MSnExp",function(object) lapply(object@spectra,mz))

##' @nord
setMethod("trimMz","Spectrum",
          function(object,mzlim,...) trimMz.Spectrum(object,mzlim))
##' @nord
setMethod("trimMz","MSnExp",
          function(object,mzlim,...) {
            object@spectra <- lapply(spectra(object),trimMz,mzlim,...)
            return(object)
          })


##' @nord
setMethod("quantify","Spectrum",
          function(object,reporters,
                   method=c("trapezoidation","max","sum")) 
          quantify.Spectrum(object,reporters,match.arg(method)))

##' @nord
setMethod("quantify","MSnExp",
          function(object,reporters,
                   method=c("trapezoidation","max","sum"),
                   verbose=TRUE) 
          quantify.MSnExp(object,reporters,match.arg(method),verbose))
         
##' @nord
setMethod("curveStats","Spectrum",
          function(object,reporters) curveStats.Spectrum(object,reporters))

##' @nord
setMethod("curveStats","MSnExp",
          function(object,reporters,verbose=TRUE) {
            ifelse(verbose,progress <- "text",progress <- "none")
            l <- llply(object@spectra,curveStats,reporters,.progress=progress)
            qdfr <- l[[1]]
            for (i in 2:length(l)) 
              qdfr <- rbind(qdfr,l[[i]])
            return(new("MSnQual",
                       qc=qdfr,
                       metadata=object@metadata,
                       process=object@process))
          })


##' @nord
setMethod("proteomicsData",
          signature(object="MSnSet"),
          function(object) object@proteomicsData)

##' @nord
setMethod("proteomicsData",
          signature(object="MSnExp"),
          function(object) object@proteomicsData)

##' @nord
setMethod("qual","MSnQual",function(object) object@qc)
##' @nord
setMethod("normalise","MSnSet",
          function(object,method=c("sum","max"))
          normalise.MSnSet(object,match.arg(method))
          )
##' @nord
setMethod("dim","MSnSet",function(x) dim(exprs(x)))
##' @nord
setMethod("dim","MSnQual",function(x) dim(qual(x)))
##' @nord
setMethod("ratios","MSnSet",function(object) ratios.MSnSet(object))

########################################
## Not sure about this...
##' @nord
setMethod("featureNames<-","MSnSet",
          function(object,value="character") object@features <- value)
##' @nord
setReplaceMethod("featureNames",
                 signature(object="MSnSet",
                           value="character"),
                 function(object, value) {
                   object@features = value
                   if (validObject(object))
                     return(object)
                 })

############################################
## Not handling Spectrum1 explicitely here, 
## as this should rarely be used by users
##' @nord
setMethod("precursorCharge<-","Spectrum2",
          function(object,value="integer") object@precursorCharge <- value)
##' @nord
setReplaceMethod("precursorCharge",
                 signature(object="Spectrum2",
                           value="integer"),
                 function(object, value) {
                   object@precursorCharge = value
                   if (validObject(object))
                     return(object)
                 })

