\name{combineFeatures}
\alias{combineFeatures}

\title{ Combines features in an 'MSnSet' object }


\description{
  This function combines the features in an
  \code{"\linkS4class{MSnSet}"} instance applying a summarisation
  function (see \code{fun} argument) to sets of features as
  defined by a factor (see \code{groupBy} argument).
  Note that the feature names are automatically updated based on the
  \code{groupBy} parameter.

  The coefficient of variations are automatically computed and collated
  to the featureData slot. See \code{cv} and \code{cv.norm} arguments
  for details.

  NB: All the functions available as \code{fun} take a \code{na.rm}
  argument. This argument is \code{FALSE} by default. This will have as
  effect that \code{NA} get propagated at the higher level. It is
  generally advised to set \code{na.rm = TRUE}. See the example below.   
}

\usage{
combineFeatures(object, groupBy, fun = c("mean", "median", "weighted.mean", "sum", "medpolish"), ..., cv = TRUE, cv.norm = "sum", verbose = TRUE)
}

\arguments{
  \item{object}{ An instance of class \code{"\linkS4class{MSnSet}"}
    whose features will be summerised. }
  \item{groupBy}{ An object of class factor defining how to summerise
    the features. }
  \item{fun}{ The summerising function. Currently, mean, median,
    weighted mean, sum and median polish are implemented, but
    user-defined functions can also be supplied. }
  \item{\dots}{ Additional arguments for the \code{fun} function. }
  \item{cv}{ A \code{logical} defining if feature coefficients of
    variation should be computed and stored as feature
    meta-data. Default is \code{TRUE}. }
  \item{cv.norm}{ A \code{character} defining how to normalise the
    feature intensitites prior to CV calculation. Default is
    \code{sum}. Use \code{none} to keep intensities as is. See
    \code{\link{featureCV}} for more details. }
  \item{verbose}{ A \code{logical} indicating whether verbose output is
    to be printed out. }
}

% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }

\value{
  A new \code{"\linkS4class{MSnSet}"} instance is returned with
  \code{ncol} (i.e. number of samples) is unchanged, but \code{nrow}
  (i.e. the number od features) is now equals to the number of levels in
  \code{groupBy}. The feature metadata (\code{featureData} slot) is
  updated accordingly and only the first occurrence of a feature in
  the original feature meta-data is kept. 
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\seealso{
  \code{\link{featureCV}}
}

\examples{
data(itraqdata)
quant <- quantify(itraqdata[11:15], method = "max", reporters = iTRAQ4)
dim(quant)
exprs(quant)
## arbitrary grouping into two groups
grp <- as.factor(c(1, 1, 2, 2, 2))
quant.comb <- combineFeatures(quant, grp, "sum")
dim(quant.comb)
exprs(quant.comb)
fvarLabels(quant.comb)

## missing data
exprs(quant)[4, 4] <- 
    exprs(quant)[2, 2] <- NA
exprs(quant)
## NAs propagate in the 115 and 117 channels 
exprs(combineFeatures(quant, grp, "sum"))
## NAs are removed before summing 
exprs(combineFeatures(quant, grp, "sum", na.rm = TRUE))
}

%% \keyword{ ~kwd1 }
