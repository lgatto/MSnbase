plot.Spectrum2 <- function(spectrum,reporters=NULL,full=FALSE) {
  mtc <- i <- NULL                             # to satisfy codetools
  xmin <- xmax <- ymin <- ymax <- fill <- NULL # to satisfy codetools
  ## plot.Spectrum2: no visible binding for global variable ‘i’
  ## ...
  if (!full & is.null(reporters))
    stop("Please provide repotrer ions if you do not want a full spectrum.")
  df <- data.frame(i=intensity(spectrum),
                   mtc=mz(spectrum))
  mainvp <- viewport(width=1,height=1,x=0.5,y=0.5)
  title <- opts(title=paste("Precursor M/Z",spectrum@precursorMz))
  p <- ggplot(df,aes(x=mtc,y=i)) + 
              labs(x="M/Z",y="Intensity (ion counts)")
  if (!is.null(reporters)) {
    rect <- data.frame(mtc=mean(df$mtc),
                       i=0,
                       xmin=reporters@mz-reporters@width,
                       xmax=reporters@mz+reporters@width,
                       ymin=0,
                       ymax=max(df$i),
                       fill=reporters@col,
                       alpha=1/3)
    p <- p + geom_rect(data=rect,
                       aes(xmin=xmin,
                           xmax=xmax,
                           ymin=ymin,
                           ymax=ymax,
                           fill=fill,
                           alpha=alpha)) 
  }
  p <- p +
    geom_point(alpha=I(1/10)) +
      geom_line() +
        opts(legend.position = "none")
  if (!is.null(reporters)) {
    if (class(reporters)!="ReporterIons")
      stop("Reporters must be of class 'ReporterIons'.")
    width <- reporters@width
    rlim1 <- min(reporters@mz)-width
    rlim2 <- max(reporters@mz)+width
    coord <- coord_cartesian(xlim=c(rlim1,rlim2))
    subvp <- viewport(width=2/3,height=1/3,x=.95,y=.92,just=c("right","top"))
    reps <- p + coord + 
      theme_gray(9) +
        labs(x = NULL, y = NULL) +
          opts(plot.margin = unit(c(1,1,0,0), "lines")) +
            scale_x_continuous(breaks=seq(rlim1,rlim2,(rlim2-rlim1)/10)) +
              opts(legend.position = "none") 
  }
  if (full) {
    print(p+title,vp=mainvp)
     if (!is.null(reporters)) 
      print(reps,vp=subvp)
    invisible(p+title)
  } else {
    print(reps+title,vp=mainvp)
    invisible(reps+title)
  }
}
