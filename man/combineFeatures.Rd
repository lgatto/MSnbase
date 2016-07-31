\name{combineFeatures}
\alias{combineFeatures}

\title{ Combines features in an \code{MSnSet} object }


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

combineFeatures(object, groupBy, fun = c("mean", "median",
"weighted.mean", "sum", "medpolish", "iPQF"), redundancy.handler =
c("unique", "multiple"), cv = TRUE, cv.norm = "sum", verbose =
isMSnbaseVerbose(), ...)

}

\arguments{

  \item{object}{ An instance of class \code{"\linkS4class{MSnSet}"}
    whose features will be summerised. }

  \item{groupBy}{ A \code{factor}, \code{character}, \code{numeric} or a
    \code{list} of the above defining how to summerise the features. The
    list must be of length \code{nrow(object)}. Each element of the list
    is a vector describing the feature mapping. If the list can be
    named, its names must match \code{fetureNames(object)}. See
    \code{redundancy.handler} for details about the latter.}

  \item{fun}{ The summerising function. Currently, mean, median,
    weighted mean, sum and median polish and iPQF (see
    \code{\link{iPQF}} for details) are implemented, but user-defined
    functions can also be supplied. }
    
  \item{redundancy.handler}{ If \code{groupBy} is a \code{list}, one of
    \code{"unique"} (default) or \code{"multiple"} (ignored otherwise)
    defining how to handle peptides that can be associated to multiple
    higher-level features (proteins) upon combination. Using
    \code{"unique"} will only consider uniquely matching features
    (features matching multiple proteins will be discarded).
    \code{"multiple"} will allow matching to multiple proteins and each
    feature will be repeatedly tallied for each possible matching
    protein. }
  
  \item{cv}{ A \code{logical} defining if feature coefficients of
    variation should be computed and stored as feature
    meta-data. Default is \code{TRUE}. }

  \item{cv.norm}{ A \code{character} defining how to normalise the
    feature intensitites prior to CV calculation. Default is
    \code{sum}. Use \code{none} to keep intensities as is. See
    \code{\link{featureCV}} for more details. }

  \item{verbose}{ A \code{logical} indicating whether verbose output is
    to be printed out. }

  \item{\dots}{ Additional arguments for the \code{fun} function. }
}


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

\references{
  iPQF: a new peptide-to-protein summarization method using peptide
  spectra characteristics to improve protein quantification. Fischer M,
  Renard BY.  Bioinformatics. 2016 Apr
  1;32(7):1040-7. doi:10.1093/bioinformatics/btv675. Epub 2015 Nov
  20. PubMed PMID:26589272.
}

\seealso{
  \code{\link{featureCV}}
}

\examples{
data(msnset)
msnset <- msnset[11:15, ]
exprs(msnset)

## arbitrary grouping into two groups
grp <- as.factor(c(1, 1, 2, 2, 2))
msnset.comb <- combineFeatures(msnset, grp, "sum")
dim(msnset.comb)
exprs(msnset.comb)
fvarLabels(msnset.comb)

## grouping with a list
grpl <- list(c("A", "B"), "A", "A", "C", c("C", "B"))
## optional naming
names(grpl) <- featureNames(msnset)
exprs(combineFeatures(msnset, grpl, fun = "sum", redundancy.handler = "unique"))
exprs(combineFeatures(msnset, grpl, fun = "sum", redundancy.handler = "multiple"))

## missing data
exprs(msnset)[4, 4] <- 
    exprs(msnset)[2, 2] <- NA
exprs(msnset)
## NAs propagate in the 115 and 117 channels 
exprs(combineFeatures(msnset, grp, "sum"))
## NAs are removed before summing 
exprs(combineFeatures(msnset, grp, "sum", na.rm = TRUE))

## using iPQF
data(msnset2)
res <- combineFeatures(msnset2,
                       groupBy = fData(msnset2)$accession,
                       redundancy.handler = "unique",
                       fun = "iPQF",
                       low.support.filter = FALSE,
                       ratio.calc = "sum",
                       method.combine = FALSE)
head(exprs(res))
}


