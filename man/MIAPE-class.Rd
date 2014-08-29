\name{MIAPE-class}
\Rdversion{1.1}
\docType{class}
\alias{MIAPE-class}
\alias{class:MIAPE}
\alias{MIAPE}

\alias{show,MIAPE-method}
\alias{msInfo,MIAPE-method}
\alias{msInfo}
\alias{expemail,MIAPE-method}
\alias{expemail}
\alias{exptitle,MIAPE-method}
\alias{exptitle}
\alias{ionSource,MIAPE-method}
\alias{ionSource}
\alias{ionSourceDetails,MIAPE-method}
\alias{ionSourceDetails}
\alias{analyser,MIAPE-method}
\alias{analyser}
\alias{analyserDetails,MIAPE-method}
\alias{analyserDetails}
\alias{analyzer,MIAPE-method}
\alias{analyzer}
\alias{analyzerDetails,MIAPE-method}
\alias{analyzerDetails}
\alias{detectorType,MIAPE-method}
\alias{detectorType}
\alias{instrumentModel,MIAPE-method}
\alias{instrumentModel}
\alias{instrumentManufacturer,MIAPE-method}
\alias{instrumentManufacturer}
\alias{instrumentCustomisations,MIAPE-method}
\alias{instrumentCustomisations}
\alias{abstract,MIAPE-method}
\alias{samples,MIAPE-method}
\alias{pubMedIds,MIAPE-method}
\alias{pubMedIds<-,MIAPE-method}
\alias{otherInfo,MIAPE-method}
\alias{expinfo,MIAPE-method}
\alias{notes,MIAPE-method}
\alias{notes<-,MIAPE-method}
\alias{combine,MIAPE,MIAPE-method}
\alias{as.MIAME.MIAPE}
\alias{coerce,MIAPE,MIAME-method}

\title{The "MIAPE" Class for Storing Proteomics Experiment Information}

\description{
  The Minimum Information About a Proteomics Experiment.
  The current implementation is based on the MIAPE-MS 2.4 document. 
}

\section{Slots}{
  \describe{
    \item{\code{title}:}{Object of class \code{character} containing a
      single-sentence experiment title. }
    \item{\code{abstract}:}{Object of class \code{character}
      containing an abstract describing the experiment. }
    \item{\code{url}:}{Object of class \code{character} containing a
      URL for the experiment. }
    \item{\code{pubMedIds}:}{Object of class \code{character} listing
      strings of PubMed identifiers of papers relevant to the dataset. }
    \item{\code{samples}:}{Object of class \code{list} containing
      information about the samples. }
    \item{\code{preprocessing}:}{Object of class \code{list} containing
      information about the pre-processing steps used on the raw data from
      this experiment. }
    \item{\code{other}:}{Object of class \code{list} containing other
      information for which none of the above slots applies. }
    \item{\code{dateStamp}:}{Object of class \code{character}, giving
      the date on which the work described was initiated; given in the
      standard 'YYYY-MM-DD' format (with hyphens). }
    \item{\code{name}:}{Object of class \code{character} containing the
      name of the (stable) primary contact person for this data set;
      this could be the experimenter, lab head, line manager, \ldots}  
    \item{\code{lab}:}{Object of class \code{character} containing the
      laboratory where the experiment was conducted. }
    \item{\code{contact}:}{Object of class \code{character} containing
      contact information for lab and/or experimenter. }
    \item{\code{email}:}{Object of class \code{character} containing
      tmail contact information for the primary contact person (see
      \code{name} above). } 
    \item{\code{instrumentModel}:}{Object of class \code{character}
      indicating the model of the mass spectrometer used to generate the
      data. }
    \item{\code{instrumentManufacturer}:}{Object of class
      \code{character} indicating the manufacturing company of the mass
      spectrometer. }
    \item{\code{instrumentCustomisations}:}{Object of class
      \code{character} describing any significant (i.e. affecting
      behaviour) deviations from the manufacturer's specification for the
      mass spectrometer. }
    \item{\code{softwareName}:}{Object of class \code{character} with
      the instrument management and data analysis package(s) name(s). }
    \item{\code{softwareVersion}:}{Object of class \code{character} with
      the instrument management and data analysis package(s) version(s). }
    \item{\code{switchingCriteria}:}{Object of class \code{character}
      describing the list of conditions that cause the switch from survey or
      zoom mode (MS1) to or tandem mode (MSn where n > 1); e.g. 'parent
      ion'' mass lists, neutral loss criteria and so on [applied for
      tandem MS only]. }
    \item{\code{isolationWidth}:}{Object of class \code{numeric}
      describing, for tandem instruments, the total width (i.e. not half
      for plus-or-minus) of the gate applied around a selected precursor ion
      m/z, provided for all levels or by MS level. }
    \item{\code{parameterFile}:}{Object of class \code{character} giving
      the location and name under which the mass spectrometer's parameter
      settings file for the run is stored, if available. Ideally this should
      be a URI+filename, or most preferably an LSID, where feasible. }
    \item{\code{ionSource}:}{Object of class \code{character} describing
      the ion source (ESI, MALDI, \ldots). }
    \item{\code{ionSourceDetails}:}{Object of class \code{character} describing
      the relevant details about the ion source. See MIAPE-MI docuement
      for more details. }
    \item{\code{analyser}:}{Object of class \code{character} describing
      the analyzer type (Quadrupole, time-of-flight, ion trap, \ldots). }
    \item{\code{analyserDetails}:}{Object of class \code{character} describing
      the relevant details about the analyzer. See MIAPE-MI document
      for more details. }
    \item{\code{collisionGas}:}{Object of class \code{character}
      describing the composition of the gas used to fragment ions in the
      collision cell. }
    \item{\code{collisionPressure}:}{Object of class \code{numeric}
      providing the pressure (in bars) of the collision gas. }
    \item{\code{collisionEnergy}:}{Object of class \code{character}
      specifying for the process of imparting a particular impetus to
      ions with a given m/z value, as they travel into the collision
      cell for fragmentation. This could be a global figure (e.g. for tandem
      TOF's), or a complex function; for example a gradient (stepped or
      continuous) of m/z values (for quads) or activation frequencies (for
      traps) with associated collision energies (given in eV). Note that
      collision energies are also provided for individual
      \code{"\linkS4class{Spectrum2}"} instances, and is the preferred
      way of accessing this data. }
    \item{\code{detectorType}:}{Object of class \code{character}
      describing the type of detector used in the machine (microchannel
      plate, channeltron, \ldots). }  
    \item{\code{detectorSensitivity}:}{Object of class \code{character}
      giving and appropriate measure of the sensitivity of the described
      detector (e.g. applied voltage). }
  }
}


\section{Methods}{
  The following methods as in \code{"\linkS4class{MIAME}"}:
  \describe{
    \item{\code{abstract(MIAPE)}:}{An accessor function for
      \code{abstract}. }
    \item{\code{expinfo(MIAPE)}:}{An accessor function for \code{name},
      \code{lab}, \code{contact}, \code{title}, and \code{url}.} 
    \item{\code{notes(MIAPE), notes(MIAPE) <- value}:}{Accessor
      functions for \code{other}. \code{notes(MIAME) <- character} 
      \emph{appends} character to notes; use \code{notes(MIAPE) <- list}
      to replace the notes entirely.} 
    \item{\code{otherInfo(MIAPE)}:}{An accessor function for
      \code{other}.} 
    \item{\code{preproc(MIAPE)}:}{An accessor function for
      \code{preprocessing}.} 
    \item{\code{pubMedIds(MIAPE), pubMedIds(MIAME) <- value}:}{Accessor
      function for \code{pubMedIds}.} 
    \item{\code{expemail(MIAPE)}:}{An accessor function for 
      \code{email} slot.}
    \item{\code{exptitle(MIAPE)}:}{An accessor function for
      \code{title} slot.}
    \item{\code{analyzer(MIAPE)}:}{ An accessor function for
      \code{analyser} slot. \code{analyser(MIAPE)} is also available.}
    \item{\code{analyzerDetails(MIAPE)}:}{ An accessor function for
      \code{analyserDetails} slot. \code{analyserDetails} is also
      available. }
    \item{\code{detectorType(MIAPE)}:}{ An accessor function for
      \code{detectorType} slot.}
    \item{\code{ionSource(MIAPE)}:}{ An accessor function for
      \code{ionSource} slot.}
    \item{\code{ionSourceDetails(MIAPE)}:}{ An accessor function for
      \code{ionSourceDetails} slot.}
    \item{\code{instrumentModel(MIAPE)}:}{ An accessor function for
      \code{instrumentModel} slot.}
    \item{\code{instrumentManufacturer(MIAPE)}:}{ An accessor function for
      \code{instrumentManufacturer} slot.}
    \item{\code{instrumentCustomisations(MIAPE)}:}{ An accessor function for
      \code{instrumentCustomisations} slot.}
    \item{\code{as(,"MIAME")}:}{Coerce the object from \code{MIAPE} to
      \code{MIAME} class. Used when converting an \code{MSnSet} into an
      \code{ExpressionSet}. }
  }
  
  MIAPE-specific methods, including MIAPE-MS meta-data:
  \describe{
    \item{\code{show(MIAPE)}:}{Displays the experiment data.}
    \item{\code{msInfo(MIAPE)}:}{Displays 'MIAPE-MS' information.}
  }
}
  
\section{Extends}{
  Class \code{"\linkS4class{MIAxE}"}, directly.
  Class \code{"\linkS4class{Versioned}"}, by class "MIAxE", distance 2.
}

\references{
  About MIAPE: \url{http://www.psidev.info/index.php?q=node/91}, and
  references therein, especially 'Guidelines for reporting the use of mass
  spectrometry in proteomics', Nature Biotechnology 26, 860-861 (2008).  
}
\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

% \note{
%   This is a temporary implementation, that will most likely be updated
%   using the AnnotedDataFrame class. MIAPE is currently just a checklist
%   of items, but may be more formally described later.
% }

\keyword{classes}
