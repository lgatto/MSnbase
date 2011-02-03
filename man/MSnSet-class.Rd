\name{MSnSet-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnSet-class}
\alias{class:MSnSet}
\alias{MSnSet}

%% TODO add methods
\alias{show,MSnSet-method}

\title{
  The "MSnSet" Class for MS Proteomics Expression Data and Meta-Data
}

\description{
  The \code{MSnSet} holds quantified expression data for MS proteomics
  data and the experimental meta-data.
  The \code{MSnSet} class is derived from the
  \code{"\linkS4class{ExpressionSet}"} class and thus benefits from many
  of the available analysis methods and requires a matrix named
  \code{exprs} as assayData member. 
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("MSnSet", exprs, ...)}.
  See also \code{"\linkS4class{ExpressionSet}"} for helpful information.
  Expression data produced from other softwares can thus make use of this
  standardized data container to benefit \code{R} and
  \code{Bioconductor} packages. Importer functions will be developed to
  stream-line the generation of \code{"MSnSet"} instances from
  third-party software. 
  
  In the frame of the \code{MSnbase} package, \code{MSnSet} instances
  are generated from \code{"\linkS4class{MSnExp}"} experiments when the
  \code{"\linkS4class{ReporterIons}"} using the \code{"quantify"} method). 
}

\section{Slots}{
  \describe{
    \item{\code{qual}:}{Object of class \code{"data.frame"} that records
      peaks data for each of the reporter ions to be used as quality
      metrics. }
    \item{\code{processingData}:}{Object of class
      \code{"\linkS4class{MSnProcess}"} that records all processing. }
    \item{\code{assayData}:}{Object of class \code{"assayData"}
      containing a matrix with equal with column number equal to
      \code{nrow(phenoData)}. \code{assayData} must contain a matrix
      \code{exprs} with rows represening features (e.g., reporters ions)
      and columns representing samples. See the \code{"\linkS4class{AssayData}"}
      class, \code{\link{exprs}} and \code{\link{assayData}} accessor
      for more details. This slot in indirectly inherited from
      \code{"\linkS4class{ExpressionSet}"}. } 
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}
      containing experimenter-supplied variables describing sample (i.e
      the individual tags for an labelled MS experiment) (indireclty 
      inherited from  \code{"\linkS4class{ExpressionSet}"}). See
      \code{\link{phenoData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }
    \item{\code{featureData}:}{Object of class
      \code{"AnnotatedDataFrame"} containing variables describing
      features (spectra in our case), e.g. identificaiton data, peptide
      sequence, identification score,... (inherited indirectly from 
      \code{"\linkS4class{ExpressionSet}"}). See
      \code{\link{featureData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }
    \item{\code{experimentData}:}{Object of class
      \code{"\linkS4class{MIAME}"}, containing details of experimental
      methods (inherited from \code{"\linkS4class{ExpressionSet}"}). See 
      \code{\link{experimentData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }
    \item{\code{proteomicsData}:}{Object of class
      \code{"\linkS4class{MIAPE}"}, containing information about the
      proteomics experiment. }
    \item{\code{annotation}:}{ not used here. }
    \item{\code{protocolData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      equipment-generated variables (inherited indirectly from
      \code{"\linkS4class{ExpressionSet}"}). See
      \code{\link{protocolData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }
    \item{\code{.__classVersion__}:}{Object of class
      \code{"\linkS4class{Versions}"} describing the versions of R,
      \code{"linkS4class{Biobase}"}, \code{"\linkS4class{eSet}"},
      \code{"\linkS4class{ExpressionSet}"} and \code{MSnSet} of the
      current instance. Intended for developer use and debugging (inherited
      from indirectly \code{"\linkS4class{ExpressionSet}"}). }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{ExpressionSet}"}, directly.
  Class \code{"\linkS4class{eSet}"}, by class "ExpressionSet", distance 2.
  Class \code{"\linkS4class{VersionedBiobase}"}, by class "ExpressionSet", distance 3.
  Class \code{"\linkS4class{Versioned}"}, by class "ExpressionSet", distance 4.
}

\section{Methods}{
  See \code{"\linkS4class{ExpressionSet}"} methods.
  %% more methods to add
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{"\linkS4class{eSet}"} and \code{"\linkS4class{ExpressionSet}"}.
}
\examples{
mzxmlfile <- dir(system.file("extdata",package="MSnbase"),
                 pattern="mzXML",full.names=TRUE)
msnexp <- readMzXMLData(mzxmlfile,verbose=FALSE)
data(iTRAQ4)
msnset <- quantify(msnexp,method="trap",reporters=iTRAQ4,verbose=FALSE)
msnset
}

\keyword{classes}
