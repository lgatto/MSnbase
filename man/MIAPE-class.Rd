\name{MIAPE-class}
\Rdversion{1.1}
\docType{class}
\alias{MIAPE-class}
\alias{class:MIAPE}
\alias{MIAPE}
\alias{proteomicsData}

\alias{show,MIAPE-method}

\title{The "MIAPE" Class for Storing Proteomics Experiment Information}

\description{
  The Minimum Information About a Proteomics Experiment. The current
  dummy implementation is composed of a free-text description of the
  experiment. The class will be updated at a later stage to be MIAPE
  compliant.
}

\section{Slots}{
  \describe{
    \item{\code{description}:}{Object of class \code{"character"}
      describing the proteomics experiment or linking to a complete
      MIAPE compliant annexe.}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      indicating the version of the \code{MIAPE} instance. Intended for
      developer use and debugging. }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  \describe{
    \item{\code{show(object)}}{ Displays the free-text description of
      the \code{MIAPE} object. }
  }
}

\references{
About MIAPE: \url{http://www.psidev.info/index.php?q=node/91} and reporting guideline: \url{http://www.psidev.info/index.php?q=node/60}
}
\author{
  Laurent Gatto <lg390@cam.ac.uk>
}
\note{
  This is a temporary implementation, that will most likely be updated
  using the AnnotedDataFrame class. MIAPE is currently just a checklist
  of items, but may be more formally described later.
}

\seealso{
  \code{"\linkS4class{MIAME}"}
}

\examples{
  new("MIAPE", description="Reference file with MIAPE meta-data.")
}

\keyword{classes}
