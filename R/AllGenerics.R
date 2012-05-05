setGeneric("clean", function(object,...) standardGeneric("clean"))
setGeneric("trimMz", function(object,mzlim,...) standardGeneric("trimMz"))
setGeneric("removePeaks", function(object,t="min",...) standardGeneric("removePeaks"))
setGeneric("removeReporters", function(object,...) standardGeneric("removeReporters"))

##setGeneric("bg.correct", function(object,bg,...) standardGeneric("bg.correct"))

setGeneric("spectra",function(object) standardGeneric("spectra"))
##setGeneric("spectra<-",function(object,value) standardGeneric("spectra<-"))
setGeneric("precursorMz",function(object) standardGeneric("precursorMz"))
setGeneric("precursorIntensity",function(object) standardGeneric("precursorIntensity"))
setGeneric("precursorCharge",function(object) standardGeneric("precursorCharge"))
setGeneric("precursorCharge<-",function(object,value) standardGeneric("precursorCharge<-"))
setGeneric("acquisitionNum",function(object) standardGeneric("acquisitionNum"))
setGeneric("precScanNum",function(object) standardGeneric("precScanNum"))
setGeneric("rtime",function(object) standardGeneric("rtime"))
setGeneric("msLevel",function(object) standardGeneric("msLevel"))
setGeneric("collisionEnergy",function(object) standardGeneric("collisionEnergy")) 
## setGeneric("peaksCount",function(object) standardGeneric("peaksCount")) ## use mzR generic
## setGeneric("header",function(object) standardGeneric("header")) ## use mzR generic
setGeneric("polarity", function(object) standardGeneric("polarity"))
setGeneric("centroided", function(object) standardGeneric("centroided"))
setGeneric("centroided<-", function(object,value) standardGeneric("centroided<-"))

setGeneric("processingData",function(object) standardGeneric("processingData"))
setGeneric("processingData<-",function(object,value) standardGeneric("processingData<-"))

## setGeneric("proteomicsData",function(object) standardGeneric("proteomicsData"))
## setGeneric("proteomicsData<-",function(object,value) standardGeneric("proteomicsData<-"))
setGeneric("msInfo",function(object) standardGeneric("msInfo"))
setGeneric("email",function(object) standardGeneric("email"))

setGeneric("mz",function(object) standardGeneric("mz"))
setGeneric("intensity",function(object) standardGeneric("intensity"))
setGeneric("tic",function(object) standardGeneric("tic"))
setGeneric("ionCount",function(object) standardGeneric("ionCount"))

setGeneric("fromFile",function(object) standardGeneric("fromFile"))

setGeneric("quantify",function(object,...) standardGeneric("quantify"))
setGeneric("curveStats",function(object,reporters,...) standardGeneric("curveStats"))
setGeneric("purityCorrect",function(object,impurities,...) standardGeneric("purityCorrect"))

setGeneric("qual",function(object) standardGeneric("qual"))
setGeneric("normalise",function(object,method,...) standardGeneric("normalise"))
setGeneric("normalize",function(object,method,...) standardGeneric("normalize"))

setGeneric("width",function(object) standardGeneric("width"))
setGeneric("reporterNames",function(object) standardGeneric("reporterNames"))
setGeneric("reporterNames<-",function(object,value) standardGeneric("reporterNames<-"))
setGeneric("reporterColours",function(object) standardGeneric("reporterColours"))
setGeneric("reporterColors",function(object) standardGeneric("reporterColors"))

setGeneric("fileNames",function(object) standardGeneric("fileNames"))
setGeneric("fileNames<-",function(object,value) standardGeneric("fileNames<-"))

setGeneric("extractPrecSpectra",function(object,prec) standardGeneric("extractPrecSpectra"))
setGeneric("extractSpectra",function(object,selected) standardGeneric("extractSpectra"))

setGeneric("multiplex",function(object) standardGeneric("multiplex"))
setGeneric("multiLabels",function(object) standardGeneric("multiLabels"))

setGeneric("plot2d",function(object,...) standardGeneric("plot2d"))
setGeneric("plotDensity",function(object,...) standardGeneric("plotDensity"))
setGeneric("plotMzDelta",function(object,...) standardGeneric("plotMzDelta"))
setGeneric("plotNA",function(object,...) standardGeneric("plotNA"))

setGeneric("writeMgfData",function(object,...) standardGeneric("writeMgfData"))

setGeneric("filterNA", function(object,...) standardGeneric("filterNA"))
setGeneric("topN", function(object,...) standardGeneric("topN"))
