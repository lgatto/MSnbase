lineplot <- function(x, ...) {
  ## x: MSnSet or matrix or data.frame
  if (class(x) == "MSnSet")
    x <- exprs(x)
  plot(0, type = "n", xlim = c(1,ncol(x)),
       ylim = range(x, na.rm = TRUE),
       xlab = "Reporter",
       ylab = "Intensity")
  grid()
  for (i in 1:nrow(x)) {
    lines(1:ncol(x), x[i,], ...)
  }
}

