\name{NAnnotatedDataFrame-class}
\Rdversion{1.1}
\docType{class}
\alias{NAnnotatedDataFrame-class}
\alias{NAnnotatedDataFrame}
\alias{class:NAnnotatedDataFrame}

\alias{multiplex,NAnnotatedDataFrame-method}
\alias{multiLabels,NAnnotatedDataFrame-method}
\alias{dim,NAnnotatedDataFrame-method}
\alias{show,NAnnotatedDataFrame-method}
\alias{multiplex}
\alias{multiLabels}


\title{Class Containing Measured Variables and Their Meta-Data
  Description for Multiplexed Experiments.
}

\description{
  An \code{NAnnotatedDataFrame} is an
  \code{"\linkS4class{AnnotatedDataFrame}"}, as defined in the 'Biobase'
  package that includes additional labels for multiplexing annotation. 
}


\section{Objects from the Class}{
  See \code{"\linkS4class{AnnotatedDataFrame}"} for object creation with
  \code{new}. Multiplexing data is defined by setting the
  \code{multiplex} and \code{multiLables} parameters.
}


\section{Slots}{
  \describe{
    \item{\code{multiplex}:}{Object of class \code{"numeric"} indicating
      the number of multiplexed samples described. }
    \item{\code{multiLabels}:}{Object of class \code{"character"}
      describing the multiplexing. }
    \item{\code{varMetadata}:}{Object of class \code{"data.frame"} with
      number of rows equal number of columns in \code{data}, and at
      least one column, named \code{labelDescription}, containing a
      textual description of each variable. Inherited from
      \code{"\linkS4class{AnnotatedDataFrame}"}. }
    \item{\code{data}:}{Object of class \code{"data.frame"}
      containing samples (rows) and measured variables
      (columns). Inherited from
      \code{"\linkS4class{AnnotatedDataFrame}"}. }
    \item{\code{dimLabels}:}{Object of class \code{"character"} of
      length 2 that provides labels for the rows and columns in the
      \code{show} method. Inherited from 
      \code{"\linkS4class{AnnotatedDataFrame}"}. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      describing the instance version. Intended for developer
      use. Inherited from  \code{"\linkS4class{AnnotatedDataFrame}"}. } 
  }
}
\section{Extends}{
Class \code{"\linkS4class{AnnotatedDataFrame}"}, directly.
Class \code{"\linkS4class{Versioned}"}, by class "AnnotatedDataFrame", distance 2.
}

\section{Methods}{
  \describe{
    \item{dim}{\code{signature(object = "NAnnotatedDataFrame")}: Returns
      the number of samples, variables and multiplex cardinality in the
      object. }
    \item{multiplex}{\code{signature(object = "NAnnotatedDataFrame")}: 
      Returns the number of multipexed samples described by the object. }
    \item{multiLabels}{\code{signature(object = "NAnnotatedDataFrame")}: 
      Returns the multiplex labels. }
    \item{show}{\code{signature(object = "NAnnotatedDataFrame")}:
      Textual description of the object. }
    }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\seealso{
  \code{"\linkS4class{AnnotatedDataFrame}"}.  
}
\examples{
df <- data.frame(x=1:3,
                 y=LETTERS[1:3],
                 row.names=paste("Sample",1:3,sep=""))
metaData <-
  data.frame(labelDescription=c(
               "Numbers",
               "Factor levels"))
mplx <- c("M1","M2")
new("NAnnotatedDataFrame",
    data=df,
    varMetadata=metaData,
    multiplex=length(mplx),
    multiLabels=mplx)
}

\keyword{classes}
