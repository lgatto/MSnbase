\name{MSmap-class}
\Rdversion{1.1}
\docType{class}
\alias{MSmap-class}
\alias{MSmap}

\alias{MSmap-method}
\alias{MSmap,OnDiskMSnExp-method}
\alias{MSmap,mzRpwiz-method}
\alias{MSmap,mzRramp-method}
\alias{MSmap,mzRraw-method}
\alias{coerce,MSmap,data.frame-method}
\alias{dim,MSmap-method}
\alias{fileName,MSmap-method}
\alias{fileNames,MSmap-method}
\alias{msLevel,MSmap-method}
\alias{msMap,MSmap-method}
\alias{msMap}
\alias{mz,MSmap-method}
\alias{mzRes}
\alias{mzRes,MSmap-method}
\alias{ncol,MSmap-method}
\alias{nrow,MSmap-method}
\alias{plot3D,MSmap-method}
\alias{plot3D}
\alias{plot,MSmap,missing-method}
\alias{rtime,MSmap-method}
\alias{show,MSmap-method}
\alias{t,MSmap-method}

\title{Class \code{MSmap}}

\description{
  A class to store mass spectrometry data maps, i.e intensities
  collected along the M/Z and retention time space during a mass
  spectrometry acquisition.
}

\section{Objects from the Class}{
  Objects can be created with the \code{MSmap} constructor. The
  constructor has the following arguments:
  
  \describe{

    \item{object}{An object created by \code{mzR::openMSfile} or an
      instance of class \code{\linkS4class{OnDiskMSnExp}}. If the latter
      contains data from multiple files, a warning will be issued and
      the first one will be used. }
      
    \item{lowMz}{A \code{numeric} of length 1 defining the lower bound
      of the M/Z range of the MS map.}

    \item{highMz}{A \code{numeric} of length 1 defining the upper bound
      of the M/Z range of the MS map. }

    \item{resMz}{The resolution along the M/Z range. }
    
    \item{hd}{An optional \code{data.frame} as produced by
      \code{mzR::header(object)}. If missing, will be computer within
      the function. Ignored when \code{object} is an
      \code{OnDiskMSnExp}.}
    
    \item{zeroIsNA}{Set 0 intensities to \code{NA}. This can be used
      to clarify the 3 dimensional plot produce by \code{plot3D}.}
  }
}


\section{Slots}{  
  \describe{
    \item{\code{call}:}{Object of class \code{"call"} - the call used to
      generate the instance. }
    \item{\code{map}:}{Object of class \code{"matrix"} containing the
      actual MS map. }
    \item{\code{mz}:}{Object of class \code{"numeric"} with the M/Z
      sampling bins. }
    \item{\code{res}:}{Object of class \code{"numeric"} storing the the M/Z
      resolution used to create the map. }
    \item{\code{rt}:}{Object of class \code{"numeric"} with the
      retention times of the map spectra. }
    \item{\code{ms}:}{Object of class \code{"numeric"} with the MS
      levels of the spectra. }    
    \item{\code{t}:}{Object of class \code{"logical"} indicating if the
    instance has been transposed. }
    \item{\code{filename}:}{Object of class \code{"character"}
      specifying the filename of the original raw MS data.  }
  }
}

\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "MSmap", to = "data.frame")}:
      convert the \code{MSmap} instance in a \code{data.frame}. Useful
      for plotting with \code{lattice} or \code{ggplot2}. }
    \item{fileName}{\code{signature(object = "MSmap")}: returns the raw
      data filename. } 
    \item{msLevel}{\code{signature(object = "MSmap")}: returns the MS
      level of the map spectra. }
    \item{msMap}{\code{signature(object = "MSmap")}: returns the actual
      map \code{matrix}. }
    \item{mz}{\code{signature(object = "MSmap", ...)}: returns the M/Z values
      of the map. Additional arguments are currently ignored. }
    \item{rtime}{\code{signature(object = "MSmap", ...)}: returns retention
      time values of the map. Additional arguments are currently ignored.}	
    \item{mzRes}{\code{signature(object = "MSmap")}: returns the
      resolution with which the sample along the M/Z range was done. }
    \item{dim}{\code{signature(x = "MSmap")}: returns the dimensions of
      the map. \code{ncol} and \code{nrow} return the number of columns
      and rows respectively. }
    \item{t}{\code{signature(x = "MSmap")}: transposes the map. }
    \item{show}{\code{signature(object = "MSmap")}: prints a summary of
      the map. }
    \item{plot}{\code{signature(x = "MSmap", allTicks = "logical")}:
      produces an image of the map using \code{lattice::levelplot}. By
      default, \code{allTicks} is \code{TRUE} and all M/Z and retention
      times ticks of drawn. If set to \code{FALSE}, only 10 ticks in
      each dimension are plotted. }
    \item{plot3D}{\code{signature(object = "MSmap", rgl = "logical")}:
      produces an three dimensional view of the map using
      \code{lattice::cloude(..., type = "h")}. If \code{rgl} is
      \code{TRUE}, the map is visualised on a \code{rgl} device and can
      be rotated with the mouse. }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{

\dontrun{
    ## downloads the data
    library("rpx")
    px1 <- PXDataset("PXD000001")
    (i <- grep("TMT.+mzML", pxfiles(px1), value = TRUE))
    mzf <- pxget(px1, i)

    ## Using an mzRpwiz object
    ## reads the data 
    ms <- openMSfile(mzf)
    hd <- header(ms)

    ## a set of spectra of interest: MS1 spectra eluted
    ## between 30 and 35 minutes retention time
    ms1 <- which(hd$msLevel == 1)
    rtsel <- hd$retentionTime[ms1] / 60 > 30 &
        hd$retentionTime[ms1] / 60 < 35

    ## the map
    M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd)

    plot(M, aspect = 1, allTicks = FALSE)
    plot3D(M)
    if (require("rgl") & interactive())
        plot3D(M, rgl = TRUE)

    ## With some MS2 spectra
    i <- ms1[which(rtsel)][1]
    j <- ms1[which(rtsel)][2]
    M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
    plot3D(M2)

    ## Using an OnDiskMSnExp object and accessors
    msn <- readMSData(mzf, mode = "onDisk")
    
    ## a set of spectra of interest: MS1 spectra eluted
    ## between 30 and 35 minutes retention time
    ms1 <- which(msLevel(msn) == 1)
    rtsel <- rtime(msn)[ms1] / 60 > 30 &
        rtime(msn)[ms1] / 60 < 35

    ## the map
    M3 <- MSmap(msn, ms1[rtsel], 521, 523, .005)
    plot(M3, aspect = 1, allTicks = FALSE)

    ## With some MS2 spectra
    i <- ms1[which(rtsel)][1]
    j <- ms1[which(rtsel)][2]
    M4 <- MSmap(msn, i:j, 100, 1000, 1)
    plot3D(M4)    
}
}

\keyword{classes}
