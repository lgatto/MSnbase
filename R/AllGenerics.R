setGeneric("clean", function(object, ...) standardGeneric("clean"))
setGeneric("estimateNoise", function(object, ...) standardGeneric("estimateNoise"))
setGeneric("pickPeaks", function(object, ...) standardGeneric("pickPeaks"))

setGeneric("trimMz", function(object, mzlim, ...) standardGeneric("trimMz"))
setGeneric("removePeaks", function(object, t="min", ...) standardGeneric("removePeaks"))
setGeneric("removeReporters", function(object, ...) standardGeneric("removeReporters"))

##setGeneric("bg.correct", function(object,bg, ...) standardGeneric("bg.correct"))

## setGeneric("peaksCount", function(object) standardGeneric("peaksCount")) ## use mzR generic
## setGeneric("header", function(object) standardGeneric("header")) ## use mzR generic

## setGeneric("proteomicsData", function(object) standardGeneric("proteomicsData"))
## setGeneric("proteomicsData<-", function(object, value) standardGeneric("proteomicsData<-"))


setGeneric("fromFile", function(object) standardGeneric("fromFile"))
setGeneric("fromFile<-", function(object, value) standardGeneric("fromFile<-"))  ## This one should remain "private"

setGeneric("curveStats", function(object,reporters, ...) standardGeneric("curveStats"))
setGeneric("purityCorrect", function(object,impurities, ...) standardGeneric("purityCorrect"))

setGeneric("qual", function(object) standardGeneric("qual"))

setGeneric("reporterNames", function(object) standardGeneric("reporterNames"))
setGeneric("reporterNames<-", function(object, value) standardGeneric("reporterNames<-"))
setGeneric("reporterColours", function(object) standardGeneric("reporterColours"))
setGeneric("reporterColors", function(object) standardGeneric("reporterColors"))

### THESE SHOULD PROBABLY BE REPLACED BY BiocGenerics::fileName?
setGeneric("fileNames", function(object, ...) standardGeneric("fileNames"))
## setGeneric("fileNames<-", function(object, value) standardGeneric("fileNames<-"))

setGeneric("extractPrecSpectra", function(object, prec) standardGeneric("extractPrecSpectra"))
setGeneric("extractSpectra", function(object, selected) standardGeneric("extractSpectra"))

setGeneric("multiplex", function(object) standardGeneric("multiplex"))
setGeneric("multiLabels", function(object) standardGeneric("multiLabels"))

setGeneric("plot2d", function(object, ...) standardGeneric("plot2d"))
setGeneric("plotDensity", function(object, ...) standardGeneric("plotDensity"))
setGeneric("plotMzDelta", function(object, ...) standardGeneric("plotMzDelta"))
setGeneric("plotNA", function(object, ...) standardGeneric("plotNA"))

setGeneric("writeMgfData", function(object, ...) standardGeneric("writeMgfData"))

setGeneric("filterZero", function(object, ...) standardGeneric("filterZero"))
setGeneric("topN", function(object, ...) standardGeneric("topN"))

## identification
setGeneric("addIdentificationData", function(object, id, ...) standardGeneric("addIdentificationData"))
setGeneric("idSummary", function(object, ...) standardGeneric("idSummary"))
setGeneric("removeNoId", function(object, ...) standardGeneric("removeNoId"))
setGeneric("removeMultipleAssignment", function(object, ...) standardGeneric("removeMultipleAssignment"))

## Feature of Interest
setGeneric("FeaturesOfInterest",
           function(fnames, description, object, ...)
           standardGeneric("FeaturesOfInterest"))
setGeneric("FoICollection", function(object, ...) standardGeneric("FoICollection"))
setGeneric("foi", function(object, ...) standardGeneric("foi"))
setGeneric("addFeaturesOfInterest", function(x, y) standardGeneric("addFeaturesOfInterest"))
setGeneric("rmFeaturesOfInterest", function(object, i) standardGeneric("rmFeaturesOfInterest"))
## setGeneric("fromIdentical", function(x, y, ...) standardGeneric("fromIdentical"))
## setGeneric("fromEqual", function(x, y, ...) standardGeneric("fromEqual"))
setGeneric("fnamesIn", function(x, y, ...) standardGeneric("fnamesIn"))


## base::trimws
.trimws.useAsDefault <- function (x, which = c("both", "left", "right"), ...)
    base::trimws(x, which, ...)

setGeneric("trimws", signature = "x",
    function(x, which = c("both", "left", "right"), ...)
        standardGeneric("trimws"),
    useAsDefault=.trimws.useAsDefault)

setGeneric("filterFile", function (object, ...) standardGeneric("filterFile"))

## isolationWindow generic is in mzR

setGeneric("bpi", function(object, ...) standardGeneric("bpi"))
setGeneric("splitByFile", function(object, f, ...) standardGeneric("splitByFile"))

## Chromatogram class:
## setGeneric("aggregationFun", function(object, ...)
##     standardGeneric("aggregationFun"))

## centroiding related
setGeneric("estimateMzResolution", function(object, ...)
    standardGeneric("estimateMzResolution"))

setGeneric("combineSpectra", function(object, ...)
           standardGeneric("combineSpectra"))

setGeneric("isolationWindowLowerMz", function(object, ...)
    standardGeneric("isolationWindowLowerMz"))
setGeneric("isolationWindowUpperMz", function(object, ...)
    standardGeneric("isolationWindowUpperMz"))
setGeneric("filterPrecursorMz", function(object, ...)
    standardGeneric("filterPrecursorMz"))
setGeneric("transformIntensity", function(object, ...)
    standardGeneric("transformIntensity"))
