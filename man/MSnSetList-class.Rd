\name{MSnSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnSetList-class}
\alias{MSnSetList}
\alias{[,MSnSetList,ANY,ANY,ANY-method}
\alias{[,MSnSetList,ANY,missing,missing-method}
\alias{[[,MSnSetList,ANY,ANY-method}
\alias{[[,MSnSetList,ANY,missing-method}
\alias{length,MSnSetList-method}
\alias{show,MSnSetList-method}
\alias{split,MSnSet,character-method}
\alias{split,MSnSet,factor-method}

\title{Class \code{"MSnSetList"}}

\description{
  A class storing lists of \code{\linkS4class{MSnSet}} instances.
}


\section{Objects from the Class}{

  Objects can be created with the

  \describe{

    \item{\code{MSnSetList(x)} constructor}{that takes a list of
      valid \code{MSnSet} instances as input \code{x}. }

    \item{\code{split(x, f, MARING)} function}{where \code{x} is a
      single \code{MSnSet}, \code{f} is either a \code{factor} or a
      \code{character} and \code{MARGIN} is either 1 or 2 when splitting
      \code{x} along the rows or the columns. When \code{f} is a
      \code{character}, it defines the feature variable (\code{MARGIN ==
      1}) or pheotype variable (\code{MARGIN == 2}) to be used for
      splitting. When \code{f} is a \code{factor}, its length will be
      matched to \code{nrow(x)} (firstly) or \code{ncol(x)} to determine
      the \code{MARGIN}.}

  }
}


\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"list"} containing valid
      \code{MSnSet} instances. }     
    \item{\code{.__classVersion__}:}{The version of the instance. For
      development purposes only. }      
  }
}
  
\section{Methods}{
  \describe{
    \item{\code{[[}}{Extracts a single \code{MSnSet} at position. }      
    \item{\code{[}}{Extracts one of more \code{MSnSets} as
      \code{MSnSetList}. }
    \item{\code{length}}{Returns the number of \code{MSnSets}. }
    \item{show}{ Display the object by printing a short summary. }
  }
}
  

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
library("pRolocdata")
data(tan2009r1)
data(tan2009r2)
## unnamed
msnl <- MSnSetList(list(tan2009r1, tan2009r2))
## named
msnl <- MSnSetList(list(A = tan2009r1, B = tan2009r2))
length(msnl)
msnl[[1]] ## an MSnSet
msnl[1]   ## an MSnSetList
}

\keyword{classes}
