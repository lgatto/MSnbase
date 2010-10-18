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
  dummy implementation is composed of a free-text description of the
  experiment in addition to the \code{"\linkS4class{MIAME}"} slots.
  The class will be updated at a later stage to be fully MIAPE
  compliant. 
}

\section{Slots}{
  \describe{
    \item{\code{proteomicsData}:}{Object of class \code{"character"}
      describing the proteomics experiment. }
    \item{\code{name}:}{Object of class \code{character} containing the
      experimenter name (inherited from \code{"linkS4class{MIAME}"}). }
    \item{\code{lab}:}{Object of class \code{character} containing the
      laboratory where the experiment was conducted  (inherited from
      \code{"linkS4class{MIAME}"}). }
    \item{\code{contact}:}{Object of class \code{character} containing
      contact information for lab and/or experimenter  (inherited from
      \code{"linkS4class{MIAME}"}). }
    \item{\code{title}:}{Object of class \code{character} containing a
      single-sentence experiment title  (inherited from
      \code{"linkS4class{MIAME}"}). }
    \item{\code{abstract}:}{Object of class \code{character} containing
      an abstract describing the experiment  (inherited from \code{"linkS4class{MIAME}"}). }
    \item{\code{url}:}{Object of class \code{character} containing a URL
      for the experiment  (inherited from \code{"linkS4class{MIAME}"}). }
    \item{\code{samples}:}{Object of class \code{list} containing
      information about the samples  (inherited from \code{"linkS4class{MIAME}"}). }
    \item{\code{hybridizations}:}{Object of class \code{list} containing
      information about the hybridizations  (inherited from \code{"linkS4class{MIAME}"}). }
    \item{\code{normControls}:}{Object of class \code{list} containing
      information about the controls such as house keeping genes
      (inherited from \code{"linkS4class{MIAME}"}). } 
    \item{\code{preprocessing}:}{Object of class \code{list} containing
      information about the pre-processing steps used on the raw data from
      this experiment (inherited from \code{"linkS4class{MIAME}"}).} 
    \item{\code{pubMedIds}:}{Object of class \code{character} listing
      strings of PubMed identifiers of papers relevant to the dataset
      (inherited from \code{"linkS4class{MIAME}"}). } 
    \item{\code{other}:}{Object of class \code{list} containing other
      information for which none of the above slots does not applies
      (inherited from \code{"linkS4class{MIAME}"}). } 
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      indicating the version of the \code{MIAPE} instance. Intended for
      developer use and debugging. }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{MIAME}"}, directly.
  Class \code{"\linkS4class{Versioned}"}, by class "MIAME", distance 2.
  Class \code{"\linkS4class{characterORMIAME}"}, by class "MIAME", distance 2.
}

\section{Methods}{
  \describe{
    \item{\code{show(object)}}{ Displays the experiment data of the
      experiment. }
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

\keyword{classes}
