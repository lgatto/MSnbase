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


setMethod("header","MSnExp",function(object) header.MSnExp(object))

setMethod("fromFile","MSnExp",function(object) unlist(eapply(assayData(object),fromFile)))

setMethod("plot",c("MSnExp","missing"),
          function(x,y,...) plot.MSnExp(x,...))

setMethod("clean","MSnExp",function(object) clean.MSnExp(object))

setMethod("[","MSnExp",
          function(x,i,j="missing",drop="missing") "[.MSnExp"(x,i))

setMethod("removePeaks","MSnExp",
          function(object,t,verbose=TRUE) removePeaks.MSnExp(object,t,verbose))

setMethod("precursorMz","MSnExp",
          function(object) {
            ## this assumes that if first spectrum
            ## has msLevel>1, all have
            if (msLevel(object)[1]>1) 
              return(unlist(eapply(assayData(object), precursorMz)))
            stop("No precursor MZ value for MS1 spectra.")
          })

setMethod("tic","MSnExp",
          function(object) sapply(assayData(object),tic))

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
              return(unlist(eapply(assayData(object), ms1scan)))
            stop("This experiment contains MS1 spectra.")
          })

setMethod("rtime","MSnExp",function(object) sapply(spectra(object),rtime))
setMethod("peaksCount","MSnExp",
          function(object) sapply(spectra(object),peaksCount))

setMethod("msLevel","MSnExp",
          function(object) unlist(eapply(assayData(object),msLevel)))

setMethod("collisionEnergy","MSnExp",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object),collisionEnergy))
            stop("No collision energy for MS1 spectra.")
          })
setMethod("intensity","MSnExp",
          function(object) eapply(assayData(object),intensity))

setMethod("mz","MSnExp",function(object) eapply(assayData(object),mz))
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
