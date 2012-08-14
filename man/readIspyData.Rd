\name{readIspyData}
\alias{readIspyData}
\concept{MSnSet}

\title{
  Reads an ispy2 result spread sheet and creates a fully featured 'MSnSet' instance.
}

\description{
  Reads an ispy2 tab-delimited spreadsheet and generates the
  corresponding \code{\linkS4class{MSnSet}} object.
}
\usage{
readIspyData(file = "ispy_results.tsv", uniquePeps = TRUE, pep = 0.05,
  na.rm = TRUE, min.int = 0, reporters = 19:23, keepAll = FALSE,
  verbose = TRUE) 
}


\arguments{
  \item{file}{A \code{character}, indicating the file name to be read
    in. Default is "ispy_results.tsv".}
  \item{uniquePeps}{A \code{logical}, indicating whether only unique
    peptides should be included. Default is TRUE. }
  \item{pep}{A \code{numeric} indicating the posterior error probability 
    threshold for peptides to be considered correctly identified. Default
    is 0.05. }
  \item{na.rm}{A \code{logical} indicating whether reporter ions
    containing one or more NA values should be excluded. Default is
    TRUE. } 
  \item{min.int}{A \code{numeric} indicating the minimal summed
    intensity threshold for reporter data to be imported. Default is
    0. Note that 'NA' values are excluded when summing the values. } 
  \item{reporters}{A \code{numeric} indicating column indices of reporter
    ions quantitation data. Default is \code{19:23} for iTRAQ 4-plex. }
  % \item{skipFillUp}{A \code{numeric} indicating which column not to fill
  %   up. Generally the column following the reporter ions quantitation
  %   data. Default is 24.}
  % \item{fillUp}{A \code{logical} specifying if empty cells should be
  %   filled up using their upper neighbour value. Default is FALSE, and
  %   warnings are thrown if NA's are found. 
  %   This can be set to TRUE for older ispy result sheets, were protein 
  %   accession numbers, description,... and other shared fields were not
  %   repeated across multiple lines. However, if one of these fields is
  %   marked NA (for instance a protein without description), the value of
  %   the previous protein will be copied (which is obviously a bug and
  %   will be fixed.)}
  \item{keepAll}{A \code{logical} that defines whether all features of
    the ispy result should be imported. If 'TRUE', 'pep', 'na.rm' and 
    'min.int' are ignored. This is equivalent to 'pep=1', 'na.rm=FALSE'
    and 'min.int=0'. Default is 'FALSE'. }
  \item{verbose}{A \code{logical} indicating whether verbose output is
    to be printed out. }
}

\value{
  An object of class \code{"\linkS4class{MSnSet}"}.
}

\references{
  Ispy is a set of perl script to analyse SILAC, 15N and MSMS data
  developed by Phil D. Charles <pdc35@cam.ac.uk> at CCP
  \url{http://www.bio.cam.ac.uk/proteomics/}.
  No ispy references published yet. 
}

\author{
  Laurent Gatto
}

\seealso{
  \code{\link{readMSData}} to import raw data.
}

\examples{
\dontrun{ispy <- readIspyData("ispy_results.tsv")}
}

\keyword{ file }
\keyword{ manip }