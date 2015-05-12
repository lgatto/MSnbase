setGeneric("width", function(x) standardGeneric("width"))

setGeneric("bin", function(object, ...) standardGeneric("bin"))
setGeneric("clean", function(object, ...) standardGeneric("clean"))
setGeneric("compareSpectra", function(object1, object2, ...) standardGeneric("compareSpectra"))
setGeneric("pickPeaks", function(object, ...) standardGeneric("pickPeaks"))
## stats::smooth already exists
setGeneric("smooth", function(x, ...) standardGeneric("smooth"))
setGeneric("trimMz", function(object, mzlim, ...) standardGeneric("trimMz"))
setGeneric("removePeaks", function(object, t="min", ...) standardGeneric("removePeaks"))
setGeneric("removeReporters", function(object, ...) standardGeneric("removeReporters"))

##setGeneric("bg.correct", function(object,bg, ...) standardGeneric("bg.correct"))

##setGeneric("spectra<-", function(object, value) standardGeneric("spectra<-"))
setGeneric("scanIndex", function(object) standardGeneric("scanIndex"))
setGeneric("precursorMz", function(object) standardGeneric("precursorMz"))
setGeneric("precursorIntensity", function(object) standardGeneric("precursorIntensity"))
setGeneric("precursorCharge", function(object) standardGeneric("precursorCharge"))
setGeneric("precursorCharge<-", function(object, value) standardGeneric("precursorCharge<-"))
setGeneric("acquisitionNum", function(object) standardGeneric("acquisitionNum"))
setGeneric("precAcquisitionNum", function(object) standardGeneric("precAcquisitionNum"))
setGeneric("precScanNum", function(object) standardGeneric("precScanNum"))
setGeneric("msLevel", function(object) standardGeneric("msLevel"))
setGeneric("collisionEnergy", function(object) standardGeneric("collisionEnergy"))
## setGeneric("peaksCount", function(object) standardGeneric("peaksCount")) ## use mzR generic
## setGeneric("header", function(object) standardGeneric("header")) ## use mzR generic
setGeneric("polarity", function(object) standardGeneric("polarity"))
setGeneric("centroided", function(object) standardGeneric("centroided"))
setGeneric("centroided<-", function(object, value) standardGeneric("centroided<-"))

setGeneric("processingData", function(object) standardGeneric("processingData"))
setGeneric("processingData<-", function(object, value) standardGeneric("processingData<-"))

## setGeneric("proteomicsData", function(object) standardGeneric("proteomicsData"))
## setGeneric("proteomicsData<-", function(object, value) standardGeneric("proteomicsData<-"))
setGeneric("msInfo", function(object) standardGeneric("msInfo"))
setGeneric("expemail", function(object) standardGeneric("expemail"))
setGeneric("exptitle", function(object) standardGeneric("exptitle"))
setGeneric("ionSource", function(object) standardGeneric("ionSource"))
setGeneric("ionSourceDetails", function(object) standardGeneric("ionSourceDetails"))
setGeneric("analyser", function(object) standardGeneric("analyser"))
setGeneric("analyzer", function(object) standardGeneric("analyzer"))
setGeneric("analyzerDetails", function(object) standardGeneric("analyzerDetails"))
setGeneric("analyserDetails", function(object) standardGeneric("analyserDetails"))
setGeneric("detectorType", function(object) standardGeneric("detectorType"))
setGeneric("instrumentModel", function(object) standardGeneric("instrumentModel"))
setGeneric("instrumentManufacturer", function(object) standardGeneric("instrumentManufacturer"))
setGeneric("instrumentCustomisations", function(object) standardGeneric("instrumentCustomisations"))

setGeneric("ionCount", function(object) standardGeneric("ionCount"))

setGeneric("fromFile", function(object) standardGeneric("fromFile"))

setGeneric("quantify", function(object, ...) standardGeneric("quantify"))
setGeneric("curveStats", function(object,reporters, ...) standardGeneric("curveStats"))
setGeneric("purityCorrect", function(object,impurities, ...) standardGeneric("purityCorrect"))

setGeneric("qual", function(object) standardGeneric("qual"))

setGeneric("reporterNames", function(object) standardGeneric("reporterNames"))
setGeneric("reporterNames<-", function(object, value) standardGeneric("reporterNames<-"))
setGeneric("reporterColours", function(object) standardGeneric("reporterColours"))
setGeneric("reporterColors", function(object) standardGeneric("reporterColors"))

### THESE SHOULD PROBABLY BE REPLACED BY BiocGenerics::fileName?
setGeneric("fileNames", function(object) standardGeneric("fileNames"))
setGeneric("fileNames<-", function(object, value) standardGeneric("fileNames<-"))

setGeneric("extractPrecSpectra", function(object, prec) standardGeneric("extractPrecSpectra"))
setGeneric("extractSpectra", function(object, selected) standardGeneric("extractSpectra"))

setGeneric("multiplex", function(object) standardGeneric("multiplex"))
setGeneric("multiLabels", function(object) standardGeneric("multiLabels"))

setGeneric("plot2d", function(object, ...) standardGeneric("plot2d"))
setGeneric("plotDensity", function(object, ...) standardGeneric("plotDensity"))
setGeneric("plotMzDelta", function(object, ...) standardGeneric("plotMzDelta"))
setGeneric("plotNA", function(object, ...) standardGeneric("plotNA"))

setGeneric("writeMgfData", function(object, ...) standardGeneric("writeMgfData"))

setGeneric("filterNA", function(object, ...) standardGeneric("filterNA"))
setGeneric("filterZero", function(object, ...) standardGeneric("filterZero"))
setGeneric("topN", function(object, ...) standardGeneric("topN"))

setGeneric("exprsToRatios", function(object, ...) standardGeneric("exprsToRatios"))
setGeneric("impute", function(object, ...) standardGeneric("impute"))

## mzR
setGeneric("chromatogram", function(object, ...) standardGeneric("chromatogram"))
setGeneric("xic", function(object, ...) standardGeneric("xic"))

## identification
setGeneric("addIdentificationData", function(object, id, ...) standardGeneric("addIdentificationData"))
setGeneric("idSummary", function(object, ...) standardGeneric("idSummary"))
setGeneric("removeNoId", function(object, ...) standardGeneric("removeNoId"))
setGeneric("removeMultipleAssignment", function(object, ...) standardGeneric("removeMultipleAssignment"))

## fragments
setGeneric("calculateFragments", function(sequence, object, ...) standardGeneric("calculateFragments"))


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


### ProtGenerics
## setGeneric("rtime", function(object, ...) standardGeneric("rtime"))
## setGeneric("tic", function(object, ...) standardGeneric("tic"))
## setGeneric("spectra", function(object, ...) standardGeneric("spectra"))
## setGeneric("mz", function(object, ...) standardGeneric("mz"))
## setGeneric("intensity", function(object, ...) standardGeneric("intensity"))
