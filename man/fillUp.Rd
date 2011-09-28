\name{fillUp}
\alias{fillUp}


\title{
  Fills up a vector
}

\description{
  This function replaces all the empty characters \code{""} and/or
  \code{NA}s with the value of the closest preceding the preceding
  non-\code{NA}/\code{""} element. The function is used to populate
  dataframe or matrice columns where only the cells of the first row in
  a set of partially identical rows are explicitly populated and the
  following are empty.
}

\usage{
fillUp(x)
}

\arguments{
  \item{x}{ a vector. }
}


\value{
  A vector as \code{x} with all empty characters \code{""} and \code{NA}
  values replaced by the preceding non-\code{NA}/\code{""} value.
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\examples{
d <- data.frame(protein=c("Prot1","","","Prot2","",""),
                peptide=c("pep11","","pep12","pep21","pep22",""),
                score=c(1:2,NA,1:3))
d
e <- apply(d,2,fillUp)
e
data.frame(e)
fillUp(d[,1])
}
