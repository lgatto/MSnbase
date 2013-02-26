fragplot <- function(x, ...) {
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


plotNA_matrix <- function(X, pNA) {
  ## no visible binding for global variable ...
  x <- value <- variable <- proteins <- y <- z <- NULL
  ## X: matrix
  ## pNA: percentrage of NAs allowed per feature
  pNA <- pNA[1]
  calcNApp <- function(x) {
    nbna <- sum(is.na(x))
    if (is.null(dim(x)))
      nbcells <- length(x)
    else 
      nbcells <- prod(dim(x))
    return(nbna/nbcells)
  }
  ocol <- apply(X,2, function(x) sum(is.na(x)))
  orow <- apply(X,1, function(x) sum(is.na(x)))
  X <- X[order(orow), order(ocol)]
  ## percentage of NA for each protein
  k <- apply(X, 1,
             function(x) 1-sum(is.na(x))/length(x))
  ## percentage of NA in data set 
  t <- sapply(1:nrow(X),
              function(x) calcNApp(X[1:x, ]))
  nkeep <- sum( k >= (1 - pNA) )
  kkeep <- 1-calcNApp(X[1:nkeep, ])  
  dfr1 <- data.frame(x = 1:length(k),
                    proteins = k,
                    data = 1 - t)
  dfr2 <- melt(dfr1, measure.vars=c("proteins", "data"))
  p <- ggplot() + 
    geom_line(data = dfr2, aes(x = x, y = value, colour = variable)) + 
      labs(x = "Protein index (ordered by data completeness)",
           y = "Data completeness") +
             theme(legend.position=c(0.23, 0.18),
                  legend.title = element_blank(),
                  legend.background = element_rect(size = 0)) +
                    scale_colour_hue(labels = c("Individual proteins", "Full dataset"),
                                     breaks = c("proteins", "data"))
  dfr0 <- data.frame(x = nrow(dfr1), y = min(dfr1$data))  
  p <- p +
    geom_point(data = dfr0, aes(x = x, y = y), alpha = 1/3) + 
      geom_text(data = dfr0, 
                aes(x = x, y = y, label = round(y, 2)),
                vjust = 1.5, size = 2.5)                 
  ## p <- p +
  ##   geom_text(data = dfr1, 
  ##             aes(x = length(proteins), y = min(data), label = round(min(data), 2)),
  ##             vjust = 1.5, size = 2.5) +
  ##               geom_point(data = dfr1, 
  ##                          aes(x = length(proteins), y = min(data)), alpha = 1/3)
  p <- p + 
    geom_text(data = data.frame(x = nkeep, y = kkeep),
              aes(x = x, y = y, label = round(y, 2)),
              hjust = -0.5, vjust = -0.5, size = 2.5) +              
                geom_point(data = data.frame(x = nkeep, y = kkeep),
                           aes(x = x, y = y), alpha = 1/3)
  
  p <- p + geom_text(data = data.frame(x = nkeep, y = (1 - pNA), z = nkeep),
                     aes(x = x, y = y, label = round(z, 2)),
                     size = 2.5, vjust = 2, hjust = 2) +
                       geom_point(data = data.frame(x = nkeep, y = (1 - pNA)),
                                  aes(x = x, y = y), alpha = 1/3)
  
  p <- p + annotate("text", label = nrow(X), x = 0, y = 1, 
                    size = 2.5, vjust = -1, alpha = 1/3)

  print(p)
  invisible(p)
}


setMethod("image", "MSnSet",
          function(x, 
                   yticks = 10,
                   x.cex.axis = .75,
                   y.cex.axis = .75,
                   ...) {
            lab <- sampleNames(x)
            x <- exprs(x)
            nc <- ncol(x)
            nr <- nrow(x)
            image(t(x), xaxt = "n", yaxt = "n", ...)
            axis(1, seq(0,1, 1/(nc - 1)),
                 labels = lab,
                 cex.axis = x.cex.axis)
            yticks <- seq(0, 1, 1/(yticks-1)) * nr
            axis(2, seq(0,1, 1/(length(yticks) - 1)),
                 labels = round(yticks, 0),
                 cex.axis = y.cex.axis)
          })
