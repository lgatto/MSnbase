\name{pSet-class}
\Rdversion{1.1}
\docType{class}
\alias{pSet-class}

\alias{dim,pSet-method}
\alias{fileNames,pSet-method}
\alias{length,pSet-method}
\alias{processingData,pSet-method}

\title{Class "pSet"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
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
Class \code{"\linkS4class{VersionedBiobase}"}, directly.
Class \code{"\linkS4class{Versioned}"}, by class "VersionedBiobase", distance 2.
}

\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "pSet")}: ... }
    \item{fileNames}{\code{signature(object = "pSet")}: ... }
    \item{length}{\code{signature(x = "pSet")}: ... }
    \item{processingData}{\code{signature(object = "pSet")}: ... }
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
showClass("pSet")
}
\keyword{classes}
