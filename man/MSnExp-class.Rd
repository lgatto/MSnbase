\name{MSnExp-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnExp-class}
\alias{acquisitionNum,MSnExp-method}
\alias{[,MSnExp-method}
\alias{clean,MSnExp-method}
\alias{collisionEnergy,MSnExp-method}
\alias{curveStats,MSnExp-method}
\alias{fromFile,MSnExp-method}
\alias{header,MSnExp-method}
\alias{intensity,MSnExp-method}
\alias{ms1scan,MSnExp-method}
\alias{msLevel,MSnExp-method}
\alias{mz,MSnExp-method}
\alias{peaksCount,MSnExp-method}
\alias{plot,MSnExp,missing-method}
\alias{polarity,MSnExp-method}
\alias{precursorCharge,MSnExp-method}
\alias{precursorMz,MSnExp-method}
\alias{proteomicsData<-,MSnExp,ANY-method}
\alias{proteomicsData<-,MSnExp,MIAPE-method}
\alias{proteomicsData,MSnExp-method}
\alias{quantify,MSnExp-method}
\alias{removePeaks,MSnExp-method}
\alias{rtime,MSnExp-method}
\alias{show,MSnExp-method}
\alias{spectra,MSnExp-method}
\alias{tic,MSnExp-method}
\alias{trimMz,MSnExp-method}

\title{Class "MSnExp"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MSnExp", versions, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"environment"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"NAnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{proteomicsData}:}{Object of class \code{"MIAPE"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{processingData}:}{Object of class \code{"MSnProcess"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{pSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "pSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "pSet", distance 3.
}

\section{Methods}{
  \describe{
    \item{acquisitionNum}{\code{signature(object = "MSnExp")}: ... }
    \item{[}{\code{signature(x = "MSnExp")}: ... }
    \item{clean}{\code{signature(object = "MSnExp")}: ... }
    \item{collisionEnergy}{\code{signature(object = "MSnExp")}: ... }
    \item{curveStats}{\code{signature(object = "MSnExp")}: ... }
    \item{fromFile}{\code{signature(object = "MSnExp")}: ... }
    \item{header}{\code{signature(object = "MSnExp")}: ... }
    \item{intensity}{\code{signature(object = "MSnExp")}: ... }
    \item{ms1scan}{\code{signature(object = "MSnExp")}: ... }
    \item{msLevel}{\code{signature(object = "MSnExp")}: ... }
    \item{mz}{\code{signature(object = "MSnExp")}: ... }
    \item{peaksCount}{\code{signature(object = "MSnExp")}: ... }
    \item{plot}{\code{signature(x = "MSnExp", y = "missing")}: ... }
    \item{polarity}{\code{signature(object = "MSnExp")}: ... }
    \item{precursorCharge}{\code{signature(object = "MSnExp")}: ... }
    \item{precursorMz}{\code{signature(object = "MSnExp")}: ... }
    \item{proteomicsData<-}{\code{signature(object = "MSnExp", value = "ANY")}: ... }
    \item{proteomicsData<-}{\code{signature(object = "MSnExp", value = "MIAPE")}: ... }
    \item{proteomicsData}{\code{signature(object = "MSnExp")}: ... }
    \item{quantify}{\code{signature(object = "MSnExp")}: ... }
    \item{removePeaks}{\code{signature(object = "MSnExp")}: ... }
    \item{rtime}{\code{signature(object = "MSnExp")}: ... }
    \item{show}{\code{signature(object = "MSnExp")}: ... }
    \item{spectra}{\code{signature(object = "MSnExp")}: ... }
    \item{tic}{\code{signature(object = "MSnExp")}: ... }
    \item{trimMz}{\code{signature(object = "MSnExp")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("MSnExp")
}
\keyword{classes}
