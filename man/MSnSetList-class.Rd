\name{MSnSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnSetList-class}
\alias{MSnSetList}
\alias{msnsets}
\alias{objlog}
\alias{[,MSnSetList,ANY,ANY,ANY-method}
\alias{[,MSnSetList,ANY,missing,missing-method}
\alias{[[,MSnSetList,ANY,ANY-method}
\alias{[[,MSnSetList,ANY,missing-method}
\alias{length,MSnSetList-method}
\alias{show,MSnSetList-method}
\alias{names,MSnSetList-method}
\alias{names<-,MSnSetList,ANY-method}
\alias{split,MSnSet,factor-method}
\alias{split,MSnSet,character-method}
\alias{unsplit,MSnSetList,factor-method}
\alias{lapply,MSnSetList-method}
\alias{sapply,MSnSetList-method}
\alias{fData,MSnSetList-method}
\alias{fData<-,MSnSetList,DataFrame-method}

\title{Storing multiple related MSnSets}

\description{
  A class for storing lists of \code{\linkS4class{MSnSet}}
  instances.
}

\details{

  There are two ways to store different sets of measurements pertaining
  an experimental unit, such as replicated measures of different
  conditions that were recorded over more than one MS
  acquisition. Without focusing on any proteomics technology in
  particular, these multiple assays can be recorded as

  \itemize{

    \item A single combined \code{MSnSet} (see the section
    \emph{Combining MSnSet instances} in the \emph{MSnbase-demo}
    section). In such cases, the different experimental (phenotypical)
    conditions are recorded as an
    \code{\linkS4class{AnnotatedDataFrame}} in the \code{phenoData}
    slots.

    Quantitative data for features that were missing in an assay are
    generally encode as missing with \code{NA} values. Alternatively,
    only features observed in all assays could be selected. See the
    \code{\link{commonFeatureNames}} functions to select only common
    features among two or more \code{MSnSet} instance.

    \item Each set of measurements is stored in an \code{MSnSet} which
    are combined into one \code{MSnSetList}. Each \code{MSnSet} elements
    can have identical or different samples and features. Unless
    compiled directly manually by the user, one would expect at least
    one of these dimensions (features/rows or samples/columns) are
    conserved (i.e. all feature or samples names are identical). See
    \code{split}/\code{unsplit} below.
    
  }
  
}

\section{Objects from the Class}{

  Objects can be created and manipluated with:

  \describe{

    \item{\code{MSnSetList(x, log, featureDAta)}}{The class constructor
      that takes a list of valid \code{MSnSet} instances as input
      \code{x}, an optional logging \code{list}, and an optional feature
      metadata \code{data.frame}.}

    \item{\code{split(x, f)}}{An \code{MSnSetList} can be created from
      an \code{\linkS4class{MSnSet}} instance. \code{x} is a single
      \code{MSnSet} and \code{f} is a \code{factor} or a
      \code{character} of length 1. In the latter case, \code{f} will be
      matched to the feature- and phenodata variable names (in that
      order). If a match is found, the respective variable is extracted,
      converted to a factor if it is not one already, and used to split
      \code{x} along the features/rows (\code{f} was a feature variable
      name) or samples/columns (\code{f} was a phenotypic variable
      name).  If \code{f} is passed as a factor, its length will be
      matched to \code{nrow(x)} or \code{ncol(x)} (in that order) to
      determine if \code{x} will be split along the features (rows) or
      sample (columns). Hence, the length of \code{f} must match exactly
      to either dimension. }

    \item{\code{unsplit(value, f)}}{The \code{unsplit} method reverses
      the effect of splitting the \code{value} \code{MSnSet} along the
      groups \code{f}.}

    \item{\code{as(x, "MSnSetList")}}{Where \code{x} is an instance of
      class \linkS4class{MzTab}. See the class documentation for
      details.}
    
  }
}

\section{Slots}{
  \describe{
    
    \item{\code{x}:}{Object of class \code{list} containing valid
      \code{MSnSet} instances. Can be extracted with the
      \code{msnsets()} accessor. }

    \item{\code{log}:}{Object of class \code{list} containing an object
      creation log, containing among other elements the call that
      generated the object. Can be accessed with \code{objlog()}. }

    \item{\code{featureData}:}{Object of class \code{DataFrame} that
    stores metadata for each object in the \code{x} slot. The number of
    rows of this \code{data.frame} must be equal to the number of items
    in the \code{x} slot and their respective (row)names must be
    identical. }
    
    \item{\code{.__classVersion__}:}{The version of the instance. For
      development purposes only. }
  }
}
  
\section{Methods}{
  \describe{
    \item{\code{"[["}}{Extracts a single \code{MSnSet} at position. }
    
    \item{\code{"["}}{Extracts one of more \code{MSnSets} as
      \code{MSnSetList}. }
    
    \item{\code{length}}{Returns the number of \code{MSnSets}. }
    
    \item{\code{names}}{Returns the names of \code{MSnSets}, if
      available. The replacement method is also available. }
    \item{\code{show}}{ Display the object by printing a short
      summary. }

    \item{\code{lapply(x, FUN, ...)}}{ Apply function \code{FUN} to each
      element of the input \code{x}. If the application of \code{FUN}
      returns and \code{MSnSet}, then the return value is an
      \code{MSnSetList}, otherwise a \code{list}}.

    \item{\code{sapply(x, FUN, ..., simplify = TRUE, USE.NAMES =
	TRUE)}}{ A \code{lapply} wrapper that simplifies the ouptut to a
	vector, matric or array is possible. See \code{?base::sapply}
	for details. }.
    
    \item{\code{fData}}{Returns the features metadata \code{featureData}
      slot. }

    \item{\code{fData<-}}{Features metadata \code{featureData}
      replacement method. }

  }
}


\seealso{The \code{\link{commonFeatureNames}} function to select common
  features among \code{MSnSet} instances. }

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
library("pRolocdata")
data(tan2009r1)
data(tan2009r2)

## The MSnSetList class
##  for an unnamed list, names are set to indices
msnl <- MSnSetList(list(tan2009r1, tan2009r2))
names(msnl)
##  a named example
msnl <- MSnSetList(list(A = tan2009r1, B = tan2009r2))
names(msnl)
msnsets(msnl)
length(msnl)
objlog(msnl)
msnl[[1]] ## an MSnSet
msnl[1]   ## an MSnSetList of length 1

## Iterating over the elements
lapply(msnl, dim) ## a list
lapply(msnl, normalise, method = "quantiles") ## an MSnSetList

fData(msnl)
fData(msnl)$X <- sapply(msnl, nrow)
fData(msnl)

## Splitting and unsplitting
##  splitting along the columns/samples
data(dunkley2006)
head(pData(dunkley2006))
(splt <- split(dunkley2006, "replicate"))
lapply(splt, dim) ## the number of rows and columns of the split elements
unsplt <- unsplit(splt, dunkley2006$replicate)
stopifnot(compareMSnSets(dunkley2006, unsplt))

##  splitting along the rows/features
head(fData(dunkley2006))
(splt <- split(dunkley2006, "markers"))
unsplt <- unsplit(splt, factor(fData(dunkley2006)$markers))
simplify2array(lapply(splt, dim))
stopifnot(compareMSnSets(dunkley2006, unsplt))

}

\keyword{classes}
