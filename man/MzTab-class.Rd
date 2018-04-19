\name{MzTab-class}
\Rdversion{1.1}
\docType{class}
\alias{MzTab-class}
\alias{class:MzTab}
\alias{MzTab}
\alias{metadata,MzTab-method}
\alias{fileName,MzTab-method}
\alias{fileNames,MzTab-method}
\alias{peptides,MzTab-method}
\alias{proteins,MzTab-method}
\alias{psms,MzTab-method}
\alias{comments}
\alias{smallMolecules}
\alias{mzTabMode}
\alias{mzTabType}
\alias{coerce,MzTab,MSnSetList-method}
\alias{show,MzTab-method}

\title{Parse \code{MzTab} files}

\description{

  The \code{MzTab} class stores the output of a basic parsing of a
  \code{mzTab} file. It contain the metadata (a \code{list}), comments
  (a \code{character} vector), and the at least of of the following data
  types: proteins, peptides, PSMs and small molecules (as
  \code{data.frames}).

  At this stage, the metadata and data are only minimally parsed. More
  specific data extraction and preparation are delegated to more
  specialised functions, such as the \code{as(., to = "MSnSetList")} and
  \code{readMzTabData} for proteomics data.

  Note that no attempts are made to verify the validitiy of the mzTab
  file. 

}

\section{Objects from the Class}{
  
  Objects can be created by calls the the constructor \code{MzTab} that
  takes a single \code{mzTab} file as input.

  The objects can subsequently be coerced to \linkS4class{MSnSetList}
  instances with \code{as(object, "MSnSetList")}. The resulting
  \code{MSnSetList} contains possibly empty \code{MSnSet} instances for
  proteins, peptide and PSMs, respectively named \code{"Proteins"},
  \code{"Peptides"} and \code{"PSMs"}.

  The assaydata slots of the two former are populated with the
  \code{protein_abundance_assay[1-n]} and
  \code{peptide_abundance_assay[1-n]} columns in the \code{mzTab}
  file. No abundance values are defined for the latter. The respective
  feature names correspond to protein accessions, peptide sequences and
  PSM identifiers, possibly made unique as by appending sequence numbers
  to duplicates.
  
}

\section{Slots}{
  \describe{
    \item{\code{Metadata}:}{Object of class \code{"list"} storing the
      metadata section. }
    \item{\code{Filename}:}{Object of class \code{"character"} storing
      the orginal file name. }
    \item{\code{Proteins}:}{Object of class \code{"data.frame"} storing
      the protein data.  }
    \item{\code{Peptides}:}{Object of class \code{"data.frame"} storing
      the peptide data. }
    \item{\code{PSMs}:}{Object of class \code{"data.frame"} storing
      the PSM data. }
    \item{\code{SmallMolecules}:}{Object of class \code{"data.frame"}
      storing the small molecules data. }
    \item{\code{Comments}:}{Object of class \code{"character"} storing
      the comments that were present in the file. }
  }
}

\section{Accessors}{  
  \describe{
    \item{metadata}{\code{signature(x = "MzTab")}: returns the meta data
      \code{list}. }
    \item{mzTabMode}{\code{signature(x = "MzTab")}: returns the mode
      (complete or summary) of the \code{mzTab} data. A shortcut for
      \code{metadata(x)$`mzTab-mode`}.}
    \item{mzTabType}{\code{signature(x = "MzTab")}: returns the type
      (quantification or identification) of the \code{mzTab} data. A
      shortcut for \code{metadata(x)$`mzTab-type`}.}    
    \item{fileName}{\code{signature(object = "MzTab")}: returns the file name
      of the original \code{mzTab} file. }
    \item{peptides}{\code{signature(object = "MzTab")}: returns the
      peptide \code{data.frame}. }
    \item{proteins}{\code{signature(object = "MzTab")}: returns the
      proteins \code{data.frame}. }
    \item{psms}{\code{signature(object = "MzTab")}: returns the
      PSMs \code{data.frame}. }
    \item{smallMolecules}{\code{signature(object = "MzTab")}: returns
      the small molecules \code{data.frame}. }
    \item{comments}{\code{signature(object = "MzTab")}: returns the
      comments. }
  }
}


\references{

  The mzTab format is a light-weight, tab-delimited file format for
  proteomics data. See https://github.com/HUPO-PSI/mzTab for details and
  specifications.

  Griss J, Jones AR, Sachsenberg T, Walzer M, Gatto L, Hartler J,
  Thallinger GG, Salek RM, Steinbeck C, Neuhauser N, Cox J, Neumann S,
  Fan J, Reisinger F, Xu QW, Del Toro N, Perez-Riverol Y, Ghali F,
  Bandeira N, Xenarios I, Kohlbacher O, Vizcaino JA, Hermjakob H. The
  mzTab data exchange format: communicating mass-spectrometry-based
  proteomics and metabolomics experimental results to a wider
  audience. Mol Cell Proteomics. 2014 Oct;13(10):2765-75. doi:
  10.1074/mcp.O113.036681. Epub 2014 Jun 30. PubMed PMID: 24980485;
  PubMed Central PMCID: PMC4189001.
  
}

\author{

  Laurent Gatto, with contributions from Richard Cotton (see
  \url{https://github.com/lgatto/MSnbase/issues/41}).

}

\examples{
## Test files from the mzTab developement repository
fls <- c("Cytidine.mzTab", "MTBLS2.mztab", 
         "PRIDE_Exp_Complete_Ac_1643.xml-mztab.txt", 
         "PRIDE_Exp_Complete_Ac_16649.xml-mztab.txt", 
         "SILAC_CQI.mzTab", "SILAC_SQ.mzTab", 
         "iTRAQ_CQI.mzTab", "iTRAQ_SQI.mzTab", 
         "labelfree_CQI.mzTab", "labelfree_SQI.mzTab", 
         "lipidomics-HFD-LD-study-PL-DG-SM.mzTab", 
         "lipidomics-HFD-LD-study-TG.mzTab")

baseUrl <- "https://raw.githubusercontent.com/HUPO-PSI/mzTab/master/examples/1_0-Proteomics-Release/"

## a list of mzTab objects
mzt <- sapply(file.path(baseUrl, fls), MzTab)
stopifnot(length(mzt) == length(fls))
mzt[[4]]

dim(proteins(mzt[[4]]))
dim(psms(mzt[[4]]))

prots4 <- proteins(mzt[[4]])
class(prots4)
prots4[1:5, 1:4]
}

\keyword{classes}
