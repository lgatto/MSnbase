##################################################################
## Methods for MSnExp class

setMethod("show",
          signature=signature(object="MSnExp"),
          function(object) {
            if (object@.cache$level > 0) {
              msnMzRange <- object@.cache$rangeMz              
              rangePrecMz <- object@.cache$rangePrecursorMz
              nPrecMz <- object@.cache$nPrecursorMz
              uPrecMz <- object@.cache$uPrecursorMz
              nrt <- object@.cache$nRtime
              rtr <- object@.cache$rangeRtime
              msLevels <- object@.cache$msLevels
              nPrecScans <- object@.cache$nPrecursorScans
              sz <- object@.cache$size
            } else {
              if (all(msLevel(object) > 1)) {
                msnPrecMz <- unname(eapply(assayData(object),precursorMz))              
                nPrecMz <- length(msnPrecMz)
                uPrecMz <- length(unique(msnPrecMz))
                rangePrecMz <- range(msnPrecMz)
                msnMzRange <- range(unname(mz(object)))
                nPrecScans <- length(unique(eapply(assayData(object),precScanNum)))
              }
              msLevels <- unique(unlist(eapply(assayData(object),msLevel)))
              sz <- sum(unlist(unname(eapply(assayData(object),object.size))))
              msnRt <- unname(eapply(assayData(object),rtime))
              nrt <- length(msnRt)
              rtr <- range(msnRt)
            }
            cat("Object of class \"",class(object),"\"\n",sep="")
            cat(" Object size in memory: ")
            if (length(assayData(object)) == 0) {
              sz <- object.size(object)
            } else {                            
              sz <- sz + object.size(object)
            }
            cat(round(sz/(1024^2),2),"Mb\n")
            cat("- - - Spectra data - - -\n")
            if (length(assayData(object)) == 0) {
              cat(" none\n")
            } else {
              cat(" MS level(s):",msLevels,"\n")
              if (all(msLevel(object) > 1)) {
                cat(" Number of MS1 acquisitions:",nPrecScans,"\n")
                cat(" Number of MSn scans:",length(ls(assayData(object))),"\n")
                cat(" Number of precursor ions:",nPrecMz,"\n")
                if (nPrecMz > 0) {
                  cat("",uPrecMz,"unique MZs\n")
                  cat(" Precursor MZ's:",paste(signif(rangePrecMz,5),collapse=" - "),"\n")
                }
                cat(" MSn M/Z range:",round(msnMzRange,2),"\n")
              } else {
                cat(" Number of MS1 scans:",length(spectra(object)),"\n")
              }
              if (nrt > 0) {
                cat(" MSn retention times:",formatRt(rtr[1]),"-",formatRt(rtr[2]),"minutes\n")
              }
            }
            show(processingData(object))
            cat("- - - Meta data  - - -\n")
            Biobase:::.showAnnotatedDataFrame(phenoData(object),
                                              labels=list(object="phenoData"))
            cat("Loaded from:\n")
            files <- processingData(object)@files
            if (length(files)>0) {
              for (i in 1:length(files)) {
                f <- basename(files[i])
                cat(" ",f,"\n")
              }
            } else {
              cat(" none\n")
            }
            Biobase:::.showAnnotatedDataFrame(protocolData(object),
                                              labels=list(object="protocolData"))
            Biobase:::.showAnnotatedDataFrame(featureData(object),
                                              labels=list(
                                                object="featureData",
                                                sampleNames="featureNames",
                                                varLabels="fvarLabels",
                                                varMetadata="fvarMetadata"))
            cat("experimentData: use 'experimentData(object)'\n")
            pmids <- pubMedIds(object)
            if (length(pmids) > 0 && all(pmids != ""))
              cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
            invisible(NULL)
          })


setMethod("plot",c("MSnExp","missing"),
          function(x,y,...) plot.MSnExp(x,...))

setMethod("plot2d",c("MSnExp"),
          function(object,z,alpha=1/3,plot=TRUE)
          plot2d.header(header(object),z,alpha,plot))

setMethod("plot2d",c("data.frame"),
          function(object,z,alpha=1/3,plot=TRUE)
          plot2d.header(object,z,alpha,plot))

setMethod("plotDensity",c("MSnExp"),
          function(object,z,log=FALSE,plot=TRUE)
          plotDensity.header(header(object),z,log,plot))

setMethod("plotDensity",c("data.frame"),
          function(object,z,log=FALSE,plot=TRUE)
          plotDensity.header(object,z,log,plot))

setMethod("plotMzDelta",c("MSnExp"),
          function(object, reporters=NULL,
                   percentage=0.1,
                   precMz=NULL,
                   precMzWidth=2,
                   bw=1,
                   xlim=c(40,200),
                   withLabels=TRUE,
                   size=2.5,
                   plot=TRUE,
                   verbose=TRUE)
          plotMzDelta.MSnExp(object, reporters, percentage,
                             precMz, precMzWidth,bw,
                             xlim, withLabels, size,
                             plot, verbose))

setMethod("clean",
          signature=signature("MSnExp"),
          function(object,verbose=TRUE) clean.MSnExp(object,verbose))

setMethod("removePeaks",signature("MSnExp"),
          function(object,t,verbose=TRUE) removePeaks.MSnExp(object,t,verbose))


setMethod("trimMz",
          signature=signature("MSnExp","numeric"),
          function(object,mzlim,...) {
            trimmed <- eapply(assayData(object),trimMz,mzlim,...)
            object@assayData <- list2env(trimmed)
            trmd <- object@processingData@trimmed
            ifelse(length(trmd)==0,
                   object@processingData@trimmed <- mzlim,
                   object@processingData@trimmed <- c(max(trmd[1],mzlim[1]),
                                                      min(trmd[2],mzlim[2])))
            object@processingData@processing <- c(object@processingData@processing,
                                                  paste("MZ trimmed [",object@processingData@trimmed[1],
                                                        "..",object@processingData@trimmed[2],"]",sep=""))
            if (object@.cache$level > 0)
              object@.cache <- setCacheEnv(assayData(object), object@.cache$level)
            if (validObject(object))
              return(object)
          })

setMethod("quantify",
          signature=signature("MSnExp"),
          function(object,
                   method = c("trapezoidation","max","sum"),
                   reporters,
                   strict = FALSE,
                   parallel = TRUE,
                   verbose = TRUE) {
            if (!inherits(reporters,"ReporterIons"))
              stop("Argument 'reporters' must inherit from 'ReporterIons' class.")
            ## this assumes that if first spectrum
            ## has msLevel>1, all have
            if (msLevel(object)[1]<2) 
              stop("No quantification for MS1 data implemented.")
            quantify.MSnExp(object, match.arg(method), reporters, strict, parallel, verbose)
          })

setMethod("curveStats","MSnExp",
          function(object,reporters,verbose=TRUE) {
            ifelse(verbose,progress <- "text",progress <- "none")
            l <- llply(object@spectra,curveStats,reporters,.progress=progress)
            qdfr <- l[[1]]
            for (i in 2:length(l)) 
              qdfr <- rbind(qdfr,l[[i]])
            return(qdfr)
          })

setMethod("extractPrecSpectra",
          signature=signature(object="MSnExp",prec="numeric"),
          function(object,prec) extractPrecSpectra.MSnExp(object,prec))

setMethod("extractSpectra",
          signature=signature(object="MSnExp",selected="logical"),
          function(object,selected) {
            msg <- c("The 'extractSpectra' function is deprecated\n",
                     "Please use the '[' subsetting operator instead.")
            .Deprecated(msg=msg)
            extractSpectra.MSnExp(object,selected)
          })

setMethod("normalise","MSnExp",
          function(object,method=c("max","sum"),...) {
            normalise.MSnExp(object,method=match.arg(method))
        })

setMethod("normalize","MSnExp",
          function(object,method="max",...)
          normalise(object,method))


setMethod("removeReporters","MSnExp",
          function(object, reporters=NULL, clean=FALSE,verbose=TRUE) {
            if (is.null(reporters))
              return(object)
            removeReporters.MSnExp(object,reporters,clean,verbose)
        })
