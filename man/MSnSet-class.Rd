\name{MSnSet-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnSet-class}
\alias{class:MSnSet}
\alias{MSnSet}

%% TODO add methods
\alias{dim,MSnSet-method}
\alias{fileNames,MSnSet-method}
%%\alias{normalise,MSnSet-method}
\alias{msInfo,MSnSet-method}
\alias{processingData,MSnSet-method}
\alias{qual,MSnSet-method}
\alias{qual}
%%\alias{ratios,MSnSet-method}
\alias{show,MSnSet-method}
\alias{purityCorrect,MSnSet-method}
\alias{purityCorrect,MSnSet,matrix-method}
\alias{meanSdPlot,MSnSet-method}
\alias{t.MSnSet}
\alias{[,MSnSet-method}
\alias{as.ExpressionSet.MSnSet}
\alias{coerce,MSnSet,ExpressionSet-method}
\alias{write.exprs}
\alias{write.exprs,MSnSet-method}
\alias{experimentData<-,MSnSet,MIAPE-method}
\alias{combine,MSnSet,MSnSet-method}
\alias{topN,MSnSet,MSnSet-method}
\alias{topN,MSnSet-method}
\alias{topN,matrix-method}
\alias{topN}
\alias{filterNA,MSnSet-method}
\alias{filterNA,matrix-method}
\alias{filterNA}
\alias{image,MSnSet-method}

%% functions
\alias{updateFvarLabels}
\alias{updateSampleNames}
\alias{updateFeatureNames}

\title{
  The "MSnSet" Class for MS Proteomics Expression Data and Meta-Data
}

\description{
  The \code{MSnSet} holds quantified expression data for MS proteomics
  data and the experimental meta-data.
  The \code{MSnSet} class is derived from the
  \code{"\linkS4class{eSet}"} class and mimics the
  \code{"\linkS4class{ExpressionSet}"} class classically used for
  microarray data.
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
      \code{"\linkS4class{eSet}"}. } 
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}
      containing experimenter-supplied variables describing sample (i.e
      the individual tags for an labelled MS experiment) (indireclty 
      inherited from  \code{"\linkS4class{eSet}"}). See
      \code{\link{phenoData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }
    \item{\code{featureData}:}{Object of class
      \code{"AnnotatedDataFrame"} containing variables describing
      features (spectra in our case), e.g. identificaiton data, peptide
      sequence, identification score,... (inherited indirectly from 
      \code{"\linkS4class{eSet}"}). See
      \code{\link{featureData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }
    \item{\code{experimentData}:}{Object of class
      \code{"\linkS4class{MIAPE}"}, containing details of experimental
      methods (inherited from \code{"\linkS4class{eSet}"}). See 
      \code{\link{experimentData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }
    \item{\code{annotation}:}{ not used here. }
    \item{\code{protocolData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      equipment-generated variables (inherited indirectly from
      \code{"\linkS4class{eSet}"}). See
      \code{\link{protocolData}} and the \code{"\linkS4class{eSet}"}
      class for more details. }
    \item{\code{.__classVersion__}:}{Object of class
      \code{"\linkS4class{Versions}"} describing the versions of R,
      the Biobase package, \code{"\linkS4class{eSet}"},
      \code{"\linkS4class{pSet}"} and \code{MSnSet} of the
      current instance. Intended for developer use and debugging (inherited
      indirectly from \code{"\linkS4class{eSet}"}). }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{eSet}"}, directly.
  Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
  Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}

\section{Methods}{
  MSnSet specific methods or over-riding it's super-class are described
  below. See also more \code{"\linkS4class{eSet}"} for
  inherited methods. 

  \describe{
    \item{dim}{\code{signature(x = "MSnSet")}: Returns the dimensions of
      object's assay data, i.e the number of samples and the number of
      features. }
    \item{fileNames}{\code{signature(object = "MSnSet")}: Access file
      names in the \code{processingData} slot. }
    \item{msInfo}{\code{signature(object = "MSnSet")}: Prints the MIAPE-MS
      meta-data stored in the \code{experimentData} slot. }
    \item{processingData}{\code{signature(object = "MSnSet")}: Access the
      \code{processingData} slot. }    
    \item{show}{\code{signature(object = "MSnSet")}: Displays object
      content as text. }
    \item{qual}{\code{signature(object = "MSnSet")}: Access the reporter
      ion peaks description. }
    \item{purityCorrect}{\code{signature(object = "MSnSet", impurities =
	"matrix")}: performs reporter ions purity correction. See
      \code{\link{purityCorrect}} documentation for more details. } 
    \item{normalise}{\code{signature(object = "MSnSet")}: Performs
      \code{MSnSet} normalisation. See \code{\link{normalise}} for more
      details. }
    \item{t}{\code{signature(x = "MSnSet")}: Returns a transposed
      \code{MSnSet} object where features are now aligned along columns 
      and samples along rows and the \code{phenoData} and
      \code{featureData} slots have been swapped. The
      \code{protocolData} slot is always dropped with a warning. }
    \item{as(,"ExpressionSet")}{\code{signature( x = "MSnSet")}: Coerce
      object from \code{MSnSet} to
      \code{\link{ExpressionSet-class}}. The \code{experimentData} slot is 
      dropped. }
    %% \item{ratios}{\code{signature(object = "MSnSet")}: ... }
    \item{write.exprs}{signature(x = "MSnSet")}{Writes expression values
      to a tab-separated file (default is \code{tmp.txt}). The
      \code{fDataCols} parameter can be used to specify which
      \code{featureData} columns (as column names, column number or
      \code{logical}) to append on the right of the expression matrix. 
      The following arguments are the same as \code{write.table}.}
    \item{combine}{signature(x = "MSnSet", y = "MSnSet", ...)}{ Combines
      2 or more \code{MSnSet} instances according to their feature names.
      Note that the \code{qual} slot and the processing information are
      silently dropped. }
    \item{topN}{signature(object = "MSnSet", groupBy, n = 3, fun, ...)}{ 
      Selects the \code{n} most intense features (typically peptides or
      spectra) out of all available for each set defined by
      \code{groupBy} (typically proteins) and creates a new instance of
      class \code{MSnSet}. If less than \code{n} features are available,
      all are selected. The \code{ncol(object)} features are summerised
      using \code{fun} (default is \code{sum}) prior to be ordered in
      decreasing order. Additional parameters can be passed to
      \code{fun} through \code{...}, for instance to control the
      behaviour of \code{topN} in case of \code{NA} values. 
      Note that the \code{qual} slot and the processing information are
      silently dropped.
      (Works also with \code{matrix} instances.)}
    \item{filterNA}{signature(object = "MSnSet", pNA = "numeric")}{ This
      methods subsets \code{object} by removing features that have more
      than \code{pNA} percent of NA values. Default \code{pNA} is 1/2.
      To remove any feature that exhibits missing data, use
      \code{pNA = 0}.
      This method also accepts \code{matrix} instances.

      See also the \code{\link{is.na.MSnSet}} and \code{\link{plotNA}}
      methods for missing data exploration.
    }
  }
}

\section{Plotting}{
  \describe{
    \item{meanSdPlot}{\code{signature(object = "MSnSet")}}{ Plots row
      standard deviations versus row means. See 
      \code{\link{meanSdPlot}} (\code{vsn} package) for more details. }    
    \item{image}{\code{signature(x = "MSnSet", yticks = "numeric",
	x.cex.axis = "numeric", y.cex.axis = "numeric", ...)}}{ Produces
      an image of expression values in the \code{x}
      object. \code{yticks} defines how many ticks should be used on the
      y axis. \code{x.cex.axis} and \code{y.cex.axis} are passed to
      calls to \code{axis} and defined the respective character
      expansion. \code{...} is passed to \code{image}.
      \item{plotNA}{\code{signature(object = "MSnSet", pNA =
	  "numeric")}}{
	Plots missing data for an \code{MSnSet} instance. \code{pNA} is a
	\code{numeric} of length 1 that specifies the percentage 
	of accepted missing data values per features. This value will be
	highlighted with a point on the figure, illustrating the overall
	percentage of NA values in the full data set and the number of
	proteins retained. Default is 1/2. See also
	\code{\link{plotNA}}.} 
    }
  }
}     

\section{Functions}{
  \describe{
    \item{updateFvarLabels}{signature(object, label, sep)}{
      This function updates \code{object}'s featureData variable labels
      by appending \code{label}. By default, \code{label} is the
      variable name and the separator \code{sep} is \code{.}.}
    \item{updateSampleNames}{signature(object, label, sep)}{
      This function updates \code{object}'s sample names by appending
      \code{label}. By default, \code{label} is the variable name and
      the separator \code{sep} is \code{.}.} 
    \item{updateFeatureNames}{signature(object, label, sep)}{
      This function updates \code{object}'s feature names by appending
      \code{label}. By default, \code{label} is the variable name and
      the separator \code{sep} is \code{.}.} 
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{"\linkS4class{eSet}"} and
  \code{"\linkS4class{ExpressionSet}"}. \code{MSnSet} quantitation
  values can be exported to a file with
  \code{\link[MSnbase]{write.exprs}}. 
}

\examples{
mzxmlfile <- dir(system.file("extdata",package="MSnbase"),
                 pattern="mzXML",full.names=TRUE)
msnexp <- readMSData(mzxmlfile,verbose=FALSE)
msnset <- quantify(msnexp,method="trap",reporters=iTRAQ4,verbose=FALSE)
msnset
}

\keyword{classes}
