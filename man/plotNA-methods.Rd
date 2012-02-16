\name{plotNA-methods}
\docType{methods}
\alias{plotNA-methods}
\alias{plotNA,MSnSet-method}
\alias{plotNA,matrix-method}
\alias{plotNA}


\title{The 'plotNA' method for 'MSnSet' missing data exploration}


\description{
  These methods produce plots that illustrate missing data. 
  The completeness of the full dataset or a set of proteins (ordered by
  increasing NA content along the x axis) is represented. 

  The methods make use the \code{ggplot2} system. An object of class
  'ggplot' is returned invisibly. 
}

\arguments{
  \item{object}{
    An object of class \code{"\linkS4class{MSnSet}"} or a \code{matrix}.
  }
  \item{pNA}{ A \code{numeric} of length 1 that specifies the percentage
    of accepted missing data values per features. This value will be
    highlighted with a point on the figure, illustrating the overall
    percentage of NA values in the full data set and the number of
    proteins retained. Default is 1/2.
  }
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnSet", ...)}}{
      Missing data for the \code{MSnSet} instance.
    }
    \item{\code{signature(object = "matrix", ...)}}{
      Missing data for the \code{matrix}.
    }
  }
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
xx <- quantify(itraqdata, reporters = iTRAQ4, verbose = FALSE)
exprs(xx)[sample(prod(dim(xx)), 120)] <- NA
plotNA(xx, pNA = 1/4)
}


\keyword{methods}
