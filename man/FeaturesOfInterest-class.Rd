\name{FeaturesOfInterest-class}
\Rdversion{1.1}
\docType{class}

\alias{FeaturesOfInterest-class}
\alias{FeaturesOfInterest}
\alias{FeaturesOfInterest-methods}
\alias{FeaturesOfInterest,character,character,missing-method}
\alias{FeaturesOfInterest,character,character,MSnSet-method}
\alias{description,FeaturesOfInterest-method}
\alias{foi,FeaturesOfInterest-method}
\alias{foi}
\alias{foi-methods}
\alias{length,FeaturesOfInterest-method}
\alias{show,FeaturesOfInterest-method}

\alias{fnamesIn}
\alias{fnamesIn-methods}
\alias{fnamesIn,FeaturesOfInterest,MSnSet-method}
\alias{fnamesIn,FeaturesOfInterest,matrix-method}

\alias{FoICollection-class}
\alias{FoICollection-methods}
\alias{FoICollection,list-method}
\alias{FoICollection,missing-method}
\alias{FoICollection}
\alias{addFeaturesOfInterest}
\alias{addFeaturesOfInterest-methods}
\alias{addFeaturesOfInterest,FeaturesOfInterest,FoICollection-method}
\alias{description,FoICollection-method}
\alias{foi,FoICollection-method}
\alias{length,FoICollection-method}
\alias{rmFeaturesOfInterest}
\alias{rmFeaturesOfInterest-methods}
\alias{rmFeaturesOfInterest,FoICollection,numeric-method}
\alias{show,FoICollection-method}

\title{Features of Interest}

\description{

  The \emph{Features of Interest} infrastructure allows to define a set
  of features of particular interest to be used/matched against existing
  data sets contained in \code{"\linkS4class{MSnSet}"}.  A specific set
  of features is stored as an \code{FeaturesOfInterest} object and a
  collection of such non-redundant instances (for example for a specific
  organism, project, ...) can be collected in a \code{FoICollection}.

}


\section{Objects from the Class}{
  
  Objects can be created with the respective \code{FeaturesOfInterest}
  and \code{FoICollection} constructors. 

}

\section{Slots}{

  \code{FeaturesOfInterest} class:

  \describe{

    \item{\code{description}:}{Object of class \code{"character"}
      describing the instance. }

    \item{\code{objpar}:}{Optional object of class \code{"list"}
      providing details about the \code{MSnSet} instance originally used
      to create the instance. See details section. }

    \item{\code{fnames}:}{Object of class \code{"character"} with the
      feature of interest names.  }

    \item{\code{date}:}{Object of class \code{"character"} with the date
      the instance was first generated. }

    \item{\code{.__classVersion__}:}{ Object of class \code{"Versions"}
      with the \code{FeaturesOfInterest} class version. Only relevant for
      development. }

  }

  \code{FoICollection} class:

  \describe{

    \item{\code{foic}:}{Object of class \code{"list"} with the
      \code{FeaturesOfInterest}. }

    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      with the \code{FoICollection} class version. Only relevant for
      development.  }
  }


}

\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}


\section{Methods}{
  \code{FeaturesOfInterest} class:
  \describe{

    \item{description}{\code{signature(object = "FeaturesOfInterest")}:
      returns the description of \code{object}. }
    
    \item{foi}{\code{signature(object = "FeaturesOfInterest")}: returns
      the features of interests. }
    
    \item{length}{\code{signature(x = "FeaturesOfInterest")}: returns
      the number of features of interest in \code{x}. }

    \item{show}{\code{signature(object = "FeaturesOfInterest")}:
      displays \code{object}. }

    \item{fnamesIn}{\code{signature(x = "FeaturesOfInterst", y =
	"MSnSet", count = "logical")}: if \code{count} is \code{FALSE}
	(default), return a logical indicating whether there is at least
	one feautre of interest present in \code{x}? Otherwise, returns
	the number of such features. Works also matrices. }
    
  }

  \code{FoICollection} class:

  \describe{

    \item{description}{\code{signature(object = "FoICollection")}:
      returns the description of \code{object}. }

    \item{foi}{\code{signature(object = "FoICollection")}: returns a
      list of \code{FeaturesOfInterest}. }

    \item{length}{\code{signature(x = "FoICollection")}: returns the
      number of \code{FeaturesOfInterest} in the collection. }

    \item{addFeaturesOfInterest}{\code{signature(x =
	"FeaturesOfInterest", y = "FoICollection")}: add the
	\code{FeaturesOfInterest} instance \code{x} to
	\code{FoICollection} \code{y}. If \code{x} is already present, a
	message is printed and \code{y} is returned unchanged. }

    \item{rmFeaturesOfInterest}{\code{signature(object =
	"FoICollection", i = "numeric")}: removes the \code{i}th
      \code{FeatureOfInterest} in the collection \code{object}. }

    \item{show}{\code{signature(object = "FoICollection")}: displays
      \code{object}. }

  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\examples{
library("pRolocdata")
data(tan2009r1)

x <- FeaturesOfInterest(description = "A traceable test set of features of interest",
                        fnames = featureNames(tan2009r1)[1:10],
                        object = tan2009r1)
x

description(x)
foi(x)

y <- FeaturesOfInterest(description = "Non-traceable features of interest",
                        fnames = featureNames(tan2009r1)[111:113])
y


xx <- FoICollection()
xx

xx <- addFeaturesOfInterest(x, xx)
xx <- addFeaturesOfInterest(y, xx)
xx

description(xx)
foi(xx)

fnamesIn(x, tan2009r1)
fnamesIn(x, tan2009r1, count = TRUE)

rmFeaturesOfInterest(xx, 1)

}

\keyword{classes}

