plot.Spectrum1 <- function(spectrum) {
  mtc <- i <- NULL # to satisfy codetools
  ## plot.Spectrum1: no visible binding for global variable ‘mtc’
  ## plot.Spectrum1: no visible binding for global variable ‘i’
  df <- data.frame(i=spectrum@intensity,
                   mtc=spectrum@mz)
  title <- opts(title=paste("Retention time",rtime(spectrum)))
  p <- ggplot(df,aes(x=mtc,y=i)) + 
              labs(x="M/Z",y="Intensity (ion counts)") +
                geom_line()
  print(p+title)
  invisible(p+title)
}
