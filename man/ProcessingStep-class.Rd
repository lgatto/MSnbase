\name{ProcessingStep-class}
\Rdversion{1.1}
\docType{class}
\concept{ProcessingStep}

\alias{ProcessingStep-class}
\alias{ProcessingStep:OnDiskMSnExp}
\alias{ProcessingStep}
\alias{show,ProcessingStep-method}

\title{Simple processing step class}

\description{
  The \code{ProcessingStep} class is a simple object to encapsule all
  relevant information of a data analysis processing step, i.e. the
  function name and all arguments.
}

\details{
  Objects of this class are mainly used to record all possible
  processing steps of an \code{\linkS4class{OnDiskMSnExp}} object for
  later \emph{lazy execution}.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("ProcessingStep",...)} or using the \code{ProcessingStep}
  constructor function.
}

\section{Slots}{
  \describe{
    \item{\code{FUN}:}{
      The function name to be executed as a character string.
    }
    \item{\code{ARGS}:}{
      A named \code{list} with all arguments to the function.
    }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\author{
  Johannes Rainer <johannes.rainer@eurac.edu>
}

\seealso{
  \code{\linkS4class{OnDiskMSnExp}}
}

\keyword{classes}
