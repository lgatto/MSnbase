\name{plotMzDelta-methods}
\docType{methods}
\alias{plotMzDelta-methods}
\alias{plotMzDelta,MSnExp-method}
\alias{plotMzDelta,mzRramp-method}
\alias{plotMzDelta}

\title{The delta m/z plot}

\description{
  The m/z delta plot illustrates the suitability of MS2 spectra for
  identification by plotting the m/z differences of the most intense
  peaks. The resulting histogram should optimally shown outstanding bars
  at amino acid residu masses. The plots have been described in Foster
  \emph{et al} 2011.

  Only a certain percentage of most intense MS2 peaks are taken into
  account to use the most significant signal. Default value is 10\% (see
  \code{percentage} argument). The difference between peaks is then
  computed for all individual spectra and their distribution is plotted
  as a histogram where single bars represent 1 m/z differences. Delta
  m/z between 40 and 200 are plotted by default, to encompass the
  residue masses of all amino acids and several common contaminants,
  although this can be changes with the \code{xlim} argument.

  In addition to the processing described above, isobaric reporter tag
  peaks (see the \code{reporters} argument) and the  precursor peak (see
  the \code{precMz} argument) can also be removed from the MS2 spectrum,
  to avoid interence with the fragment peaks.

  Note that figures in Foster \emph{et al} 2011 have been produced and
  optimised for centroided data. Application of the plot as is for data
  in profile mode has not been tested thoroughly, although the example
  below suggest that it might work.
  
  The methods make use the \code{ggplot2} system. An object of class
  \code{ggplot} is returned invisibly.

  Most of the code for plotMzDelta has kindly been contributed by
  Guangchuang Yu. 
}

\arguments{
  \item{object}{ An object of class \code{MSnExp} or \code{mzRramp}
    (from the \code{mzR} package) containing MS2 spectra. } 
  \item{reporters}{ An object of class class
    \code{"\linkS4class{ReporterIons}"} that defines which reporter ion
    peaks to set to 0. The default value \code{NULL} leaves the spectra
    as they are. }
  \item{subset}{ A numeric between 0 and 1 to use a subset of 
    \code{object}'s MS2 spectra. } 
  \item{percentage}{ The percentage of most intense peaks to be used for
    the plot. Default is 0.1.}
  \item{precMz}{ A \code{numeric} of length one or \code{NULL}
    default. In the latter (and preferred) case, the precursor m/z
    values are extracted from the individual MS2 spectra using the
    \code{\link{precursorMz}} method. }
  \item{precMzWidth}{ A \code{numeric} of length 1 that specifies the 
    width around the precursor m/z where peaks are set to 0. Default is
    2. }
  \item{bw}{ A \code{numeric} specifying the bandwith to be used to bin
    the delta m/z value to plot the histogram. Default if 1. See
    \code{\link{geom_histogram}} for more details. }
  \item{xlim}{ A \code{numeric} of length 2 specifying the range of
    delta m/z to plot on the histogram. Default is \code{c(40,200)}. }
  \item{withLabels}{ A \code{logical} defining if amino acid residue
    labels are plotted on the figure. Default is \code{TRUE}. }
  \item{size}{ A \code{numeric} of length 1 specifying the font size of
    amino acids lables. Default is 2.5.}
  \item{plot}{ A \code{logical} of length 1 that defines whether the
    figure should be plotted on the active device. Default is
    \code{TRUE}. Note that the \code{ggplot} object is always returned
    invisibly. }
  \item{verbose}{ A \code{logical} of length 1 specifying whether
    textual output and a progress bar illustration the progress of data
    processing should be printed. Default is \code{TRUE} }
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", ...)}}{
      Plots and (invisibly) returns the m/z delta histogram.
    }
  }
}

\seealso{
  The \code{\link{plotDensity}} and \code{\link{plot2d}} methods for
  other QC plots. 
}

\references{
  Foster JM, Degroeve S, Gatto L, Visser M, Wang R, Griss J, Apweiler R,
  Martens L. "A posteriori quality control for the curation and reuse of
  public proteomics data." \emph{Proteomics}, 2011
  Jun;11(11):2182-94. doi:10.1002/pmic.201000602. Epub 2011 May 2.
  PMID: 21538885 
}

\author{
  Laurent Gatto <lg390@cam.ac.uk> and Guangchuang Yu
}

\examples{
mzdplot <- plotMzDelta(itraqdata,
                       subset = 0.5,
                       reporters = iTRAQ4,
                       verbose = FALSE, plot = FALSE)
## let's retrieve peptide sequence information
## and get a table of amino acids
peps <- as.character(fData(itraqdata)$PeptideSequence)
aas <- unlist(strsplit(peps,""))
## table of aas
table(aas)
## mzDelta plot
print(mzdplot)
}

\keyword{methods}
