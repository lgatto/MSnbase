\name{NAnnotatedDataFrame-class}
\Rdversion{1.1}
\docType{class}
\alias{NAnnotatedDataFrame-class}
\alias{NAnnotatedDataFrame}
\alias{class:NAnnotatedDataFrame}

\title{Class "NAnnotatedDataFrame"}

\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NAnnotatedDataFrame", data, varMetadata, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{multiplex}:}{Object of class \code{"numeric"} indicating
      the number of multiplexed samples described. }
    \item{\code{multiLabels}:}{Object of class \code{"character"}
      describing the multiplexing. }
    \item{\code{varMetadata}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{dimLabels}:}{Object of class \code{"character"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{AnnotatedDataFrame}"}, directly.
Class \code{"\linkS4class{Versioned}"}, by class "AnnotatedDataFrame", distance 2.
}
\section{Methods}{
No methods defined with class "NAnnotatedDataFrame" in the signature.
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
showClass("NAnnotatedDataFrame")
}
\keyword{classes}
