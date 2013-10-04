plot_Spectrum1 <- function(spectrum,
                           centroided.=FALSE,
                           plot=TRUE,
                           w1) {
  mtc <- mz(spectrum)
  i <- intensity(spectrum)
  if (centroided.) {
    if (missing(w1))
      w1 <- max(mtc)/500
    dfr <- data.frame(i=i,mtc=mtc,width=w1)      
    p <- ggplot(dfr,aes(x=mtc,y=i,width=width)) + 
      geom_bar(stat="identity",position="identity")
  } else {
    dfr <- data.frame(i=i,mtc=mtc)
    p <- ggplot(dfr,aes(x=mtc,y=i)) + geom_line()
  }
  title <- theme(title=paste("Retention time",rtime(spectrum)))
  p <- p + labs(x="M/Z",y="Intensity") + title
  if (plot)
    print(p+title)
  invisible(p+title)
}
