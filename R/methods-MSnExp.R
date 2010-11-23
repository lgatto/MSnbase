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
          function(object,t,verbose=TRUE) removePeaks.MSnExp(object,t,verbose))
setMethod("bg.correct","MSnExp",
          function(object,bg,verbose=TRUE)
          bg.correct.MSnExp(object,bg=-1,verbose))
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
setMethod("ms1scan","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), ms1scan))
            stop("This experiment contains MS1 spectra.")
          })
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
