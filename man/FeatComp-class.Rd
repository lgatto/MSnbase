\name{FeatComp-class}
\Rdversion{1.1}
\docType{class}
\alias{FeatComp-class}
\alias{compfnames-methods}
\alias{compfnames,MSnSet,MSnSet-method}
\alias{compfnames,list,missing-method}
\alias{compfnames}
\alias{show,FeatComp-method}
\alias{names,FeatComp-method}
\alias{common,FeatComp-method}
\alias{common,methods}
\alias{common}
\alias{unique1,FeatComp-method}
\alias{unique1,methods}
\alias{unique1}
\alias{unique2,FeatComp-method}
\alias{unique2,methods}
\alias{unique2}

\title{Class \code{"FeatComp"}}

\description{
  Comparing feature names of two comparable \code{MSnSet} instances.
}

\section{Objects from the Class}{

  Objects can be created with \code{compfnames}. The method compares the
  feature names of two objects of class \code{"MSnSet"}. It prints a
  summary matrix of common and unique feature names and invisibly
  returns a list of \code{FeatComp} instances.

  The function will compute the common and unique features for all
  feature names of the two input objects (\code{featureNames(x)} and
  \code{feautreNames(y)}) as well as distinct subsets as defined in the
  \code{fcol1} and \code{fcol2} feautre variables.

}

\section{Slots}{
  \describe{

    \item{\code{name}:}{Object of class \code{"character"} defining the
      name of the compared features. By convention, \code{"all"} is used
      when all feature names are used; otherwise, the respective levels of
      the feature variables \code{fcol1} and \code{fcol2}. }

    \item{\code{common}:}{Object of class \code{"character"} with the common
      feature names. }

    \item{\code{unique1}:}{Object of class \code{"character"} with the
      features unique to the first \code{MSnSet} (\code{x} in
      \code{compfname}). }

    \item{\code{unique2}:}{Object of class \code{"character"} with the
      features unique to the seconn \code{MSnSet} (\code{y} in
      \code{compfname}). }

    \item{\code{all}:}{Object of class \code{"logical"} defining if all
      features of only a subset were compared. One expects that
      \code{name == "all"} when \code{all} is \code{TRUE}. }
    
  }
}

\section{Methods}{

  Accessors \code{names}, \code{common}, \code{unique1} and
  \code{unique2} can be used to access the respective \code{FeatComp}
  slots.

  \describe{

    \item{compfnames}{\code{signature(x = "MSnSet", y = "MSnSet", fcol1
	= "character", fcol2 = "character", simplify = "logical",
	verbose = "logical")}: creates the \code{FeatComp} comparison
      object for instances \code{x} and \code{y}. The feature
      variables to be considered to details feature comparison can be
      defined by \code{fcol1} (default is \code{"markers"} and
      \code{fcol2} for \code{x} and \code{y} respectively). Setting
      either to \code{NULL} will only consider all feature names; in
      such case, of \code{simplify} is \code{TRUE} (default), an
      \code{FeatComp} object is returned instead of a list of length
      1. The \code{verbose} logical controls if a summary table needs
      to be printed (default is \code{TRUE}). }

    \item{compfnames}{\code{signature(x = "list", y = "missing", ...)}:
      when \code{x} is a list of \code{MSnSet} instances,
      \code{compfnames} is applied to all element pairs of
      \code{x}. Additional parameters \code{fcol1}, \code{fcol2},
      \code{simplify} and \code{verbose} are passed to the pairwise
      comparison method. }
    
    \item{show}{\code{signature(object = "FeatComp")}: prints a summary
      of the \code{object}. }

  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk> and Thomas Naake
}

\seealso{
  \code{\link{averageMSnSet}} to compuate an average \code{MSnSet}.
}

\examples{
library("pRolocdata")
data(tan2009r1)
data(tan2009r2)
x <- compfnames(tan2009r1, tan2009r2)
x[[1]]
x[2:3]
head(common(x[[1]]))

data(tan2009r3)
tanl <- list(tan2009r1, tan2009r2, tan2009r3)
xx <- compfnames(tanl, fcol1 = NULL)
length(xx)
tail(xx)

all.equal(xx[[15]],
          compfnames(tan2009r2, tan2009r3, fcol1 = NULL))
str(sapply(xx, common))
}

\keyword{classes}
