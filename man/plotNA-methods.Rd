\name{plotNA-methods}
\docType{methods}
\alias{plotNA-methods}
\alias{plotNA,MSnSet-method}
\alias{plotNA,matrix-method}
\alias{plotNA}
\alias{is.na.MSnSet}


\title{Exploring missing data in 'MSnSet' instances}

\description{
  These methods produce plots that illustrate missing data.

  \code{is.na} returns the expression matrix of it \code{MSnSet}
  argument as a matrix of logicals referring whether the corresponding
  cells are \code{NA} or not. It is generally used in conjunction with
  \code{table} and \code{image} (see example below).

  The \code{plotNA} method produces plots that illustrate missing data.
  The completeness of the full dataset or a set of proteins (ordered by
  increasing NA content along the x axis) is represented.
  The methods make use the \code{ggplot2} system. An object of class
  'ggplot' is returned invisibly.
}


\section{Methods}{
  \describe{
    \item{is.na}{\code{signature(x = "MSnSet")}
      Returns the a matrix of logicals of dimensions \code{dim(x)}
      specifiying if respective values are missing in the
      \code{MSnSet}'s expression matrix.
    }
    \item{plotNA}{\code{signature(object = "MSnSet", pNA = "numeric")}
      Plots missing data for an \code{MSnSet} instance. \code{pNA} is a
      \code{numeric} of length 1 that specifies the percentage
      of accepted missing data values per features. This value will be
      highlighted with a point on the figure, illustrating the overall
      percentage of NA values in the full data set and the number of
      proteins retained. Default is 1/2.
    }
  }
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  See also the \code{\link{filterNA}} method to filter out features with
  a specified proportion if missing values.
}

\examples{
data(msnset)
exprs(msnset)[sample(prod(dim(msnset)), 120)] <- NA

head(is.na(msnset))
table(is.na(msnset))
image(msnset)

plotNA(msnset, pNA = 1/4)
}


\keyword{methods}
