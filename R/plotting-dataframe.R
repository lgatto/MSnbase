plot2d.header <- function(object, ## MSnExp header
                          z=c("tic","file","peaks.count","charge"),
                          alpha, plot) {
  z <- match.arg(z)
  stopifnot(c("retention.time","precursor.mz",z) %in% names(object))
  peaks.count <- charge <- retention.time <- precursor.mz <- NULL # to satisfy codetools
  p <- ggplot(object,aes(retention.time,precursor.mz)) + labs(colour=z)
  switch(z,
         tic = p <- p + geom_point(aes(colour=tic), alpha=alpha),
         peaks.count = p <- p + geom_point(aes(colour=peaks.count), alpha=alpha),
         file = p <- p + geom_point(aes(colour=as.factor(file)), alpha=alpha),
         charge = p <- p + geom_point(aes(colour=as.factor(charge)), alpha=alpha))
  if (plot)
    print(p)
  invisible(p)
}

plotDensity.header <- function(object, ## MSnExp header
                               z=c("precursor.mz","peaks.count","tic"),
                               log, plot) {
  z <- match.arg(z)
  stopifnot(c("charge",z) %in% names(object))
  peaks.count <- charge <- tic <- precursor.mz <- NULL # to satisfy codetools
  object$charge <- as.factor(object$charge)
  switch(z,
         precursor.mz = p <- ggplot(object,aes(precursor.mz)) + xlab("Precursor M/Z"),
         peaks.count = p <- ggplot(object,aes(peaks.count)) + xlab("Peaks count"),
         tic = p <- ggplot(object,aes(tic)) + xlab("Total ion current"))
  p <- p + geom_freqpoly(aes(group=charge,colour=charge))
  if (log)
    p <- p + scale_x_log10()
  if (plot)
    print(p)
  invisible(p)
}
