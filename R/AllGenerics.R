setGeneric("clean", function(object) standardGeneric("clean"))
setGeneric("trimMz", function(object,mzlim,...) standardGeneric("trimMz"))
setGeneric("removePeaks", function(object,t="min",...) standardGeneric("removePeaks"))
##setGeneric("bg.correct", function(object,bg,...) standardGeneric("bg.correct"))

setGeneric("spectra",function(object) standardGeneric("spectra"))
##setGeneric("spectra<-",function(object,value) standardGeneric("spectra<-"))
setGeneric("precursorMz",function(object) standardGeneric("precursorMz"))
setGeneric("precursorCharge",function(object) standardGeneric("precursorCharge"))
setGeneric("precursorCharge<-",function(object,value) standardGeneric("precursorCharge<-"))
setGeneric("acquisitionNum",function(object) standardGeneric("acquisitionNum"))
setGeneric("ms1scan",function(object) standardGeneric("ms1scan"))
setGeneric("rtime",function(object) standardGeneric("rtime"))
setGeneric("peaksCount",function(object) standardGeneric("peaksCount"))
setGeneric("msLevel",function(object) standardGeneric("msLevel"))
setGeneric("collisionEnergy",function(object) standardGeneric("collisionEnergy"))
setGeneric("header",function(object) standardGeneric("header"))

setGeneric("processingData",function(object) standardGeneric("processingData"))
setGeneric("processingData<-",function(object,value) standardGeneric("processingData<-"))

setGeneric("proteomicsData",function(object) standardGeneric("proteomicsData"))
setGeneric("proteomicsData<-",function(object,value) standardGeneric("proteomicsData<-"))

setGeneric("mz",function(object) standardGeneric("mz"))
setGeneric("intensity",function(object) standardGeneric("intensity"))
setGeneric("tic",function(object) standardGeneric("tic"))

setGeneric("fromFile",function(object) standardGeneric("fromFile"))

setGeneric("quantify",function(object,reporters,method,...) standardGeneric("quantify"))
setGeneric("curveStats",function(object,reporters,...) standardGeneric("curveStats"))

setGeneric("qual",function(object) standardGeneric("qual"))
setGeneric("normalise",function(object,method) standardGeneric("normalise"))
setGeneric("ratios",function(object) standardGeneric("ratios"))

setGeneric("reporterNames",function(object) standardGeneric("reporterNames"))

setGeneric("fileNames",function(object) standardGeneric("fileNames"))
setGeneric("fileNames<-",function(object,value) standardGeneric("fileNames<-"))

