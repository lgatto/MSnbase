plot.Spectrum2 <- function(spectrum,reporters=NULL,full=FALSE) {
  if (!full & is.null(reporters))
    stop("Please provide repotrer ions if you do not want a full spectrum.")
  df <- data.frame(i=spectrum@intensity,
                   mtc=spectrum@mz)
  mainvp <- viewport(width=1,height=1,x=0.5,y=0.5)
  title <- opts(title=paste("Precursor M/Z",spectrum@precursorMz))
  p <- ggplot(df,aes(x=mtc,y=i)) + 
              labs(x="M/Z",y="Intensity (ion counts)") +
              geom_line()
  if (!is.null(reporters)) {
    if (class(reporters)!="ReporterIons")
      stop("Reporters must be of class \"ReporterIons\".")
    width <- reporters@width
    rlim1 <- min(reporters@mz)-width
    rlim2 <- max(reporters@mz)+width
    coord <- coord_cartesian(xlim=c(rlim1,rlim2))
    subvp <- viewport(width=2/3,height=1/3,x=.95,y=.92,just=c("right","top"))
    reps <- p + coord + geom_point(alpha=I(1/10)) +
      theme_gray(9) +
        labs(x = NULL, y = NULL) +
          opts(plot.margin = unit(c(1,1,0,0), "lines")) +
            scale_x_continuous(breaks=seq(rlim1,rlim2,(rlim2-rlim1)/10))
  }
  if (full) {
    print(p+title,vp=mainvp)
    invisible(p+title)
    if (!is.null(reporters)) 
      print(reps,vp=subvp)
  } else {
    print(reps+title,vp=mainvp)
    invisible(reps+title)
  }
}
