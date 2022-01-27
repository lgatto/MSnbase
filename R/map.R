.MSmap <- setClass("MSmap",
                   slots = c(
                       call = "call",
                       map = "matrix",
                       mz = "numeric",
                       res = "numeric",
                       rt = "numeric",
                       ms = "numeric",
                       t = "logical",
                       filename = "character"))

setGeneric("MSmap", function(object, ...) standardGeneric("MSmap"))
setGeneric("msMap", function(object, ...) standardGeneric("msMap"))
setGeneric("mzRes", function(object, ...) standardGeneric("mzRes"))
setGeneric("plot3D", function(object, ...) standardGeneric("plot3D"))

setMethod("dim", "MSmap", function(x) dim(x@map))
setMethod("nrow", "MSmap", function(x) nrow(x@map))
setMethod("ncol", "MSmap", function(x) ncol(x@map))
setMethod("fileName", "MSmap", function(object, ...) object@filename)
setMethod("fileNames", "MSmap", function(object, ...) object@filename)
setMethod("msMap", "MSmap", function(object) object@map)
setMethod("rtime", "MSmap", function(object) object@rt)
setMethod("msLevel", "MSmap", function(object) object@ms)
setMethod("mz", "MSmap", function(object) object@mz)
setMethod("mzRes", "MSmap", function(object) object@res)

setMethod("t", "MSmap", function(x) {
    .MSmap(call = x@call,
           map = t(x@map),
           mz = x@mz,
           res = x@res,
           rt = x@rt,
           t = !x@t,
           filename = x@filename)
})


setMethod("show", "MSmap",
          function(object) {
              cat("Object of class \"", class(object), "\"\n", sep = "")
              cat(" Map [", nrow(object), ", ", ncol(object), "]\n", sep = "")
              msg <- c("rt" = paste("  Retention time:", formatRt(min(object@rt)),
                           "-", formatRt(max(object@rt)), "\n"),
                       "mz" = paste("  M/Z: ", min(object@mz), " - ",
                           max(object@mz), " (res ", mzRes(object),
                           ")\n", sep = ""))
              if (object@t) {
                  cat("  [1]", msg["mz"], sep = "")
                  cat("  [2]", msg["rt"], sep = "")
              } else {
                  cat("  [1]", msg["rt"], sep = "")
                  cat("  [2]", msg["mz"], sep = "")
              }
          })

as.data.frame.MSmap <- function(x) as(x, "data.frame")

setAs("MSmap", "data.frame",
      function(from) {
          .int <- as.numeric(msMap(from))
          .rt <- rtime(from)/60
          .mz <- mz(from)
          .ms <- msLevel(from)
          data.frame(intensity = as.numeric(.int),
                     rt = rep(.rt, length(.mz)),
                     mz = rep(.mz, each = length(.rt)),
                     ms = rep(.ms, length(.mz)))
      })

setMethod("plot", c("MSmap", "missing"),
          function(x, y, allTicks = TRUE, ...){
              m <- msMap(x)
              colnames(m) <- mz(x)
              rownames(m) <- formatRt(rtime(x))
              scales.set <- list()
              if (!allTicks) {
                  i <- as.integer(quantile(1:length(rtime(x)),
                                           seq(0, 1, length.out = 10)))
                  j <- as.integer(quantile(1:length(mz(x)),
                                           seq(0, 1, length.out = 10)))
                  scales.set <-
                      list(x = list(at = i, labels = formatRt(rtime(x)[i])),
                           y = list(at = j, labels = mz(x)[j]))
              }
              levelplot(log10(m),
                        scales = scales.set,
                        xlab = "Retention time",
                        ylab = "M/Z", ...)
          })

setMethod("plot3D", "MSmap",
          function(object, rgl = FALSE) {
              dd <- as(object, "data.frame")
              if (rgl) {
                  if (!require("rgl"))
                      stop("The 'rgl' package needed. Install it with 'install.packages(\"rgl\")'.")
                  rgl::plot3d(dd$mz, dd$rt, dd$intensity, , type = "h",
                         xlab = "M/Z", ylab = "Retention time", zlab = "")
              } else {
                  ms <- NULL ## get rid of 'no visible global function definition' note
                  par.set <- list(box.3d = list(lwd=.2))
                  cloud(intensity ~ mz + rt , data = dd,
                        type="h",
                        scales= list(
                            arrows=FALSE,
                            cex=.65,
                            draw=TRUE),
                        aspect=c(.8, 1),
                        group = ms,
                        zoom = 1,
                        par.settings = par.set,
                        axis.line = list(col = "transparent"),
                        xlab="M/Z", ylab="Retention time", zlab=NULL)
              }
          })


setClassUnion("mzRraw", c("mzRpwiz"))

setMethod("MSmap", "mzRraw",
          function(object, scans, lowMz, highMz, resMz, hd,
                   zeroIsNA = FALSE) {
              if (missing(hd))
                  hd <- header(object)
              ms1 <- which(hd$msLevel == 1)
              if (missing(scans))
                  scans <- ms1
              .call <- match.call()
              map <- mzR::get3Dmap(object, scans, lowMz, highMz, resMz)
              if (zeroIsNA)
                  map[map == 0] <- NA
              mz <- seq(lowMz, highMz, resMz)
              rt <- hd$retentionTime[scans]
              .MSmap(call = .call, map = map,
                     mz = mz, res = resMz,
                     rt = rt, ms = hd$msLevel[scans],
                     t = FALSE,
                     filename = mzR::fileName(object))
          })


setMethod("MSmap", "OnDiskMSnExp",
          function(object, scans, lowMz, highMz, resMz, hd = NULL,
                   zeroIsNA = FALSE) {
              ms1 <- which(msLevel(object) == 1)
              if (missing(scans))
                  scans <- ms1
              .call <- match.call()
              fn <- fileNames(object)
              if (length(fn) > 1)
                  warning("Using first file to build the map.")
              fh <- mzR::openMSfile(fn)
              map <- mzR::get3Dmap(fh, scans, lowMz, highMz, resMz)
              if (zeroIsNA)
                  map[map == 0] <- NA
              mz <- seq(lowMz, highMz, resMz)
              rt <- rtime(object)[scans]
              .MSmap(call = .call, map = map,
                     mz = mz, res = resMz,
                     rt = rt, ms = msLevel(object)[scans],
                     t = FALSE,
                     filename = fn)
          })
