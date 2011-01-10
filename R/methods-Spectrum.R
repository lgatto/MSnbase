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
setMethod("ms1scan","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@ms1scan)
            stop("This is already an MS1 spectrum.")
          })
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
setMethod("tic","Spectrum",function(object) sum(object@intensity))

setMethod("trimMz","Spectrum",
          function(object,mzlim,...) trimMz.Spectrum(object,mzlim))
setMethod("quantify","Spectrum",
          function(object,reporters,
                   method=c("trapezoidation","max","sum")) 
          quantify.Spectrum(object,reporters,match.arg(method)))
setMethod("curveStats","Spectrum",
          function(object,reporters) curveStats.Spectrum(object,reporters))

setMethod("precursorCharge<-","Spectrum",
          function(object,value="integer") object@precursorCharge <- value)
setReplaceMethod("precursorCharge",
                 signature(object="Spectrum",
                           value="integer"),
                 function(object, value) {
                   object@precursorCharge = value
                   if (validObject(object))
                     return(object)
                 })

setMethod("fromFile","Spectrum", function(object) object@fromFile)
