\name{MIAPE-class}
\Rdversion{1.1}
\docType{class}
\alias{MIAPE-class}
\alias{class:MIAPE}
\alias{MIAPE}

\alias{show,MIAPE-method}

\title{The "MIAPE" Class for Storing Proteomics Experiment Information}

\description{
  The Minimum Information About a Proteomics Experiment. The current
  dummy implementation is composed of a free-text description slot.
  Additional experiment met data can be set/found in
  \code{"experimentData"} slot (see\code{"\linkS4class{MIAME}"} class
  for more details).
  The class will be updated at a later stage to be fully MIAPE
  compliant. 
}

\section{Slots}{
  \describe{
    \item{\code{description}:}{Object of class \code{"character"}
      describing the proteomics experiment. }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{MIAxE}"}, directly.
  Class \code{"\linkS4class{Versioned}"}, by class "MIAxE", distance 2.
}

\section{Methods}{
  \describe{
    \item{\code{show(object)}}{ Displays the experiment data of the
      experiment. }
  }
}

\examples{
my.miape <- new("MIAPE",
  description="Relevant data about my proteomics experiment comes here.") 
my.miape
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

\keyword{classes}
