plot.MSnExp <- function(object,reporters,full=FALSE) {
  if (any(msLevel(object)==1)) {
    stop("MSnExp object contains MS1 level spectra.\n  Plotting only supported for MSn spectra")
  }
  spectra <- spectra(spectra)
  ints <- unlist(sapply(spectra, function(x) x@intensity))
  mzs <- unlist(sapply(spectra, function(x) x@mz))
  l <- unlist(sapply(spectra, function(x) length(x@mz)))
  n <- rep(1:length(l),l)
  dfr <- data.frame(i=ints,mz=mzs,n=n)
  pmz <- paste(unique(unlist(sapply(spectra,function(x) x@precursorMz))),collapse=",")
  title <- opts(title=paste("Precursor M/Z",pmz))
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
  print(p)
  invisible(p)
}

