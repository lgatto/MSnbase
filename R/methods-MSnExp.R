##################################################################
## Methods for MSnExp class
setMethod("show",
          signature=signature(object="MSnExp"),
          function(object) {
            cat("Object of class \"",class(object),"\"\n",sep="")
            cat(" Object size in memory: ")
            if (length(assayData(object))==0) {
              sz <- 0
            } else {
              sz <- sum(sapply(assayData(object),object.size)) + object.size(object)
            }
            cat(round(sz/(1024^2),2),"Mb\n")
            cat("- - - Spectra data - - -\n")
            if (sz==0) {
              cat(" none\n")
            } else {
              msnLevels <- unique(msLevel(object))
              cat(" MSn level(s):",msnLevels,"\n")
              if (all(msLevel(object)>1)) {
                cat(" Number of MS1 acquisitions:",length(unique(ms1scan(object))),"\n")
                cat(" Number of MS2 scans:",length(ls(assayData(object))),"\n")
                msnPrecMz <- precursorMz(object)
                nbPrecIons <- length(msnPrecMz)
                cat(" Number of precursor ions:",nbPrecIons,"\n")
                if (nbPrecIons>0) {
                  cat("",length(unique(msnPrecMz)),"unique MZs\n")
                  cat(" Precursor MZ's:",paste(signif(range(msnPrecMz),5),collapse=" - "),"\n")
                }
                msnMzRange <- round(range(mz(object)),2)
                cat(" MSn M/Z range:",msnMzRange,"\n")
              } else {
                cat(" Number of MS1 scans:",length(spectra(object)),"\n")
              }
              msnRt <- rtime(object)
              if (length(msnRt)>0) {
                rtr <- range(msnRt)
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
            cat("proteomicsData: use 'proteomicsData(object)'\n")
            pmids <- pubMedIds(object)
            if (length(pmids) > 0 && all(pmids != ""))
              cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
          })


setMethod("plot",c("MSnExp","missing"),
          function(x,y,...) plot.MSnExp(x,...))

setMethod("clean",
          signature=signature("MSnExp"),
          function(object,verbose=TRUE) clean.MSnExp(object,verbose))

setMethod("removePeaks","MSnExp",
          function(object,t,verbose=TRUE) removePeaks.MSnExp(object,t,verbose))


setMethod("ms1scan","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(unlist(eapply(assayData(object), ms1scan)))
            stop("This experiment contains MS1 spectra.")
          })

setMethod("trimMz",
          signature=signature("MSnExp","numeric"),
          function(object,mzlim,...) {
            trimmed <- eapply(assayData(object),trimMz,mzlim,...)
            object@assayData <- list2env(trimmed)
            return(object)
          })

setMethod("quantify",
          signature=signature("MSnExp","character"),
          function(object,
                   method=c("trapezoidation","max","sum"),
                   reporters,
                   verbose=TRUE) {
            if (!inherits(reporters,"ReporterIons"))
              stop("Argument 'reporters' must inherit from 'ReporterIons' class.")
            ## this assumes that if first spectrum
            ## has msLevel>1, all have
            if (msLevel(object)[1]<2) 
              stop("No quantification for MS1 data implemented.")
            quantify.MSnExp(object,reporters,match.arg(method),verbose)
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
          function(object,selected) extractSpectra.MSnExp(object,selected))
