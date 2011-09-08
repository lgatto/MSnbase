##################################################################
## Methods for Spectrum class and children
setMethod("initialize",
          "Spectrum",
          function(.Object,...,mz,intensity,peaksCount) {
            if (!missing(mz) & !missing(intensity)) {
              .Object <- callNextMethod(.Object,
                                        ...,
                                        intensity=intensity,
                                        mz=mz,
                                        peaksCount=length(mz))
            } else if (!missing(mz) | !missing(intensity)) {
              stop("'mz' and 'intensity' or none required.")              
            } else {
              .Object <- callNextMethod(.Object,...)
            }
            if (validObject(.Object))
              .Object
          })

setMethod("show","Spectrum",
          function(object) {
            if (msLevel(object)==1) show.Spectrum1(object)
            else show.Spectrum2(object)
            invisible(NULL)
          })

setMethod("plot",c("Spectrum","missing"),
          function(x,y,...) {
            if (msLevel(x)==1) plot.Spectrum1(x,...)
            else plot.Spectrum2(x,...)
          })

setMethod("clean",
          signature=signature("Spectrum"),
          function(object) clean.Spectrum(object))

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

setMethod("precursorIntensity","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@precursorIntensity)
            stop("No precursor data for MS1 spectra.")
          })

setMethod("acquisitionNum","Spectrum",function(object) object@acquisitionNum)
setMethod("ms1scan","Spectrum",
          function(object) {
            if (msLevel(object)>1) 
              return(object@ms1scan)
            stop("This is already an MS1 spectrum.")
          })
setMethod("rtime","Spectrum",function(object) object@rt)

setMethod("peaksCount",
          signature("Spectrum","missing"),
          function(object) object@peaksCount)
## no singature("Spectrum","numeric") -- does not make sense

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

setMethod("trimMz",
          signature=signature("Spectrum","numeric"),
          function(object,mzlim,...) trimMz.Spectrum(object,mzlim))

setMethod("quantify",
          signature=signature("Spectrum"),
          function(object,
                   method=c("trapezoidation","max","sum"),
                   reporters,
                   strict=FALSE) {
            if (!inherits(reporters,"ReporterIons"))
              stop("Argument 'reporters' must inherit from 'ReporterIons' class.")
            quantify.Spectrum(object,match.arg(method),reporters,strict)
          })

setMethod("curveStats","Spectrum",
          function(object,reporters) curveStats.Spectrum(object,reporters))

setReplaceMethod("precursorCharge",
                 signature(object="Spectrum",
                           value="integer"),
                 function(object, value) {
                   object@precursorCharge <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("fromFile","Spectrum", function(object) object@fromFile)

setMethod("polarity","Spectrum",
          function(object) {
            if (msLevel(object)==1) 
              return(object@polarity)
            stop("No polarity for MS2 spectra.")
          })

setAs("Spectrum", "data.frame",
      function (from)
      data.frame(mz=mz(from),
                 i=intensity(from))
      )

as.data.frame.Spectrum <- function(x, row.names=NULL, optional=FALSE, ...)
  as(x, "data.frame")

setMethod("centroided","Spectrum",function(object) object@centroided)

setReplaceMethod("centroided",
                 signature(object="Spectrum",
                           value="logical"),
                 function(object, value) {
                   object@centroided <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("normalise","Spectrum",
          function(object,method=c("max","sum"),...) {
            normalise.Spectrum(object,method=match.arg(method))
        })

setMethod("normalize","Spectrum",
          function(object,method=c("max","sum"),...)
          normalise(object,method=method))

setMethod("removeReporters","Spectrum",
          function(object,reporters=NULL,clean=FALSE) {
            if (msLevel(object)>1) 
              return(removeReporters.Spectrum2(object,reporters,clean))
            stop("No reporters to remove for MS1 spectra.")            
          })
