plot.MSnExp <- function(object,reporters,full=FALSE,plot=TRUE) {
  i <- NULL # to satisfy codetools
  ## plot.MSnExp: no visible binding for global variable ‘i’
  spectraList <- spectra(object)
  ints <- unlist(sapply(spectraList, function(x) x@intensity))
  mzs <- unlist(sapply(spectraList, function(x) x@mz))
  l <- unlist(sapply(spectraList, function(x) length(x@mz)))
  n <- rep(1:length(l),l)
  dfr <- data.frame(i=ints,mz=mzs,n=n)
  colnames(dfr) <- c("i","mz","n")
  if (all(msLevel(object)>1)) {
    pmz <- paste(unique(unlist(sapply(spectraList,function(x) round(precursorMz(x),2)))),collapse=",")
    title <- opts(title=paste("Precursor M/Z",pmz))
  } else {
    rtm <- paste(formatRt(range(rtime(object))),collapse=" - ")
    title <- opts(title=paste("Retention time",rtm))
    full <- TRUE
  }
  p <- ggplot(data=dfr,aes(x=mz,y=i)) +
    geom_line()+
      facet_grid(n~.) +
##        geom_point(alpha=I(1/10)) +
      title
  if (!full) {
    if (class(reporters)!="ReporterIons")
      stop("Reporters must be of class \"ReporterIons\".")
    width <- reporters@width
    rlim1 <- min(reporters@mz)-width
    rlim2 <- max(reporters@mz)+width
    reps <- coord_cartesian(xlim=c(rlim1,rlim2)) 
    breaks <- scale_x_continuous(breaks=seq(rlim1,rlim2,(rlim2-rlim1)/10))
    p <- p + reps + breaks
##    + geom_vline(xintercept=c(reporters@mz+reporters@width,
##                   reporters@mz-reporters@width),col="grey")
  }
  if (plot)
    print(p)
  invisible(p)
}

