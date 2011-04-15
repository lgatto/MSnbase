plot.Spectrum1 <- function(spectrum,
                           centroided=FALSE,
                           w1) {
  if (missing(w1))
    w1 <- max(mtc)/500
  mtc <- mz(spectrum)
  i <- intensity(spectrum)
  dfr <- data.frame(i=i,mtc=mtc,width=w1)
  title <- opts(title=paste("Retention time",rtime(spectrum)))
  if (centroided) {
    p <- ggplot(df,aes(x=mtc,y=i,width=width)) + 
      geom_bar(stat="identity",position="identity")
  } else {
    p <- ggplot(df,aes(x=mtc,y=i)) + geom_line()
  }
  p <- p + labs(x="M/Z",y="Intensity (ion counts)")  
  print(p+title)
  invisible(p+title)
}
