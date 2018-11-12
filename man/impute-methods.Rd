\name{impute-methods}
\docType{methods}
\alias{impute}
\alias{impute-methods}
\alias{impute,MSnSet-method}
\alias{imputeMethods}
\alias{naset}


\title{ Quantitative proteomics data imputation }

\description{ The \code{impute} method performs data imputation on an
  \code{MSnSet} instance using a variety of methods (see below).  The
  imputation and the parameters are logged into the
  \code{processingData(object)} slot.

  Users should proceed with care when imputing data and take precautions
  to assure that the imputation produce valid results, in particular
  with naive imputations such as replacing missing values with 0.
}


\details{

  There are two types of mechanisms resulting in missing values in
  LC/MSMS experiments.

  \itemize{
    \item{Missing values resulting from absence of detection of a
      feature, despite ions being present at detectable concentrations.
      For example in the case of ion suppression or as a result from the
      stochastic, data-dependent nature of the MS acquisition
      method. These missing value are expected to be randomly
      distributed in the data and are defined as missing at random (MAR)
      or missing completely at random (MCAR). }

    \item{Biologically relevant missing values resulting from the
      absence of the low abundance of ions (below the limit of detection
      of the instrument). These missing values are not expected to be
      randomly distributed in the data and are defined as missing not at
      random (MNAR). }

  }


  MNAR features should ideally be imputed with a left-censor method,
  such as \code{QRILC} below. Conversely, it is recommended to use host
  deck methods such nearest neighbours, Bayesian missing value
  imputation or maximum likelihood methods when values are missing at
  random.

  Currently, the following imputation methods are available:

  \describe{

    \item{MLE}{Maximum likelihood-based imputation method using the EM
      algorithm. Implemented in the \code{norm::imp.norm} function. See
      \code{\link[norm]{imp.norm}} for details and additional
      parameters. Note that here, \code{...} are passed to the
      \code{\link[norm]{em.norm}} function, rather to the actual
      imputation function \code{imp.norm}. }

    \item{bpca}{Bayesian missing value imputation are available, as
    implemented in the and \code{pcaMethods::pca} functions. See
    \code{\link[pcaMethods]{pca}} for details and additional
    parameters. }

    \item{knn}{Nearest neighbour averaging, as implemented in the
      \code{impute::impute.knn} function. See
      \code{\link[impute]{impute.knn}} for details and additional
      parameters. }

    \item{QRILC}{A missing data imputation method that performs the
      imputation of left-censored missing data using random draws from a
      truncated distribution with parameters estimated using quantile
      regression. Implemented in the \code{imputeLCMD::impute.QRILC}
      function. See \code{\link[imputeLCMD]{impute.QRILC}} for details
      and additional parameters. }

    \item{MinDet}{Performs the imputation of left-censored missing data
      using a deterministic minimal value approach. Considering a
      expression data with \emph{n} samples and \emph{p} features, for
      each sample, the missing entries are replaced with a minimal value
      observed in that sample. The minimal value observed is estimated as
      being the q-th quantile (default \code{q = 0.01}) of the observed
      values in that sample. Implemented in the
      \code{imputeLCMD::impute.MinDet} function. See
      \code{\link[imputeLCMD]{impute.MinDet}} for details and additional
      parameters. }

    \item{MinProb}{Performs the imputation of left-censored missing data
      by random draws from a Gaussian distribution centred to a minimal
      value. Considering an expression data matrix with \emph{n} samples
      and \emph{p} features, for each sample, the mean value of the
      Gaussian distribution is set to a minimal observed value in that
      sample. The minimal value observed is estimated as being the q-th
      quantile (default \code{q = 0.01}) of the observed values in that
      sample. The standard deviation is estimated as the median of the
      feature standard deviations. Note that when estimating the
      standard deviation of the Gaussian distribution, only the
      peptides/proteins which present more than 50\% recorded values are
      considered. Implemented in the \code{imputeLCMD::impute.MinProb}
      function. See \code{\link[imputeLCMD]{impute.MinProb}} for details
      and additional parameters. }

    \item{min}{Replaces the missing values by the smallest non-missing
      value in the data. }

    \item{zero}{Replaces the missing values by 0.}

    \item{mixed}{A mixed imputation applying two methods (to be defined
      by the user as \code{mar} for values missing at random and
      \code{mnar} for values missing not at random, see example) on two
      M[C]AR/MNAR subsets of the data (as defined by the user by a
      \code{randna} logical, of length equal to \code{nrow(object)}). }

    \item{nbavg}{Average neighbour imputation for fractions collected
      along a fractionation/separation gradient, such as sub-cellular
      fractions. The method assumes that the fraction are ordered along
      the gradient and is invalid otherwise.

      Continuous sets \code{NA} value at the beginning and the end of
      the quantitation vectors are set to the lowest observed value in
      the data or to a user defined value passed as argument \code{k}.
      Them, when a missing value is flanked by two non-missing
      neighbouring values, it is imputed by the mean of its direct
      neighbours. A stretch of 2 or more missing values will not be
      imputed. See the example below.

    }

    \item{none}{No imputation is performed and the missing values are
      left untouched. Implemented in case one wants to only impute value
      missing at random or not at random with the \code{mixed} method.}

  }

  The \code{naset} \code{\linkS4class{MSnSet}} is an real quantitative
  data where quantitative values have been replaced by \code{NA}s. See
  \code{script/naset.R} for details.

}



\section{Methods}{
  \describe{

    \item{\code{signature(object = "MSnSet", method, ...)}}{ This method
      performs data imputation on the \code{object} \code{MSnSet}
      instance using the \code{method} algorithm. \code{...} is used to
      pass parameters to the imputation function. See the respective
      methods for details and additional parameters.}

  }
}

\examples{
data(naset)
## table of missing values along the rows
table(fData(naset)$nNA)
## table of missing values along the columns
pData(naset)$nNA

## non-random missing values
notna <- which(!fData(naset)$randna)
length(notna)
notna

impute(naset, method = "min")

if (require("imputeLCMD")) {
    impute(naset, method = "QRILC")
    impute(naset, method = "MinDet")
}

if (require("norm"))
    impute(naset, method = "MLE")

impute(naset, "mixed",
       randna = fData(naset)$randna,
       mar = "knn", mnar = "QRILC")

## neighbour averaging

x <- naset[1:4, 1:6]
exprs(x)[1, 1] <- NA ## min value
exprs(x)[2, 3] <- NA ## average
exprs(x)[3, 1:2] <- NA ## min value and average
## 4th row: no imputation
exprs(x)

exprs(impute(x, "nbavg"))
}

\references{

Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor
Hastie, Robert Tibshirani, David Botstein and Russ B.  Altman, Missing
value estimation methods for DNA microarrays Bioinformatics (2001) 17
(6): 520-525.

Oba et al., A Bayesian missing value estimation method for gene
expression profile data, Bioinformatics (2003) 19 (16): 2088-2096.

Cosmin Lazar (2015). imputeLCMD: A collection of methods for
left-censored missing data imputation. R package version
2.0. \url{http://CRAN.R-project.org/package=imputeLCMD}.

Lazar C, Gatto L, Ferro M, Bruley C, Burger T. Accounting for the
Multiple Natures of Missing Values in Label-Free Quantitative Proteomics
Data Sets to Compare Imputation Strategies. J Proteome Res. 2016 Apr
1;15(4):1116-25. doi: 10.1021/acs.jproteome.5b00981. PubMed PMID:
26906401.

}

\author{
  Laurent Gatto and Samuel Wieczorek
}

\keyword{methods}
