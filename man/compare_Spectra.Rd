\name{compare_Spectra}
\alias{compare_Spectra}
\title{calculate similarity between spectra (between their intensity profile)}
\usage{
compare_Spectra(x, y, fun = c("common", "cor", "dotproduct"), ...)
}
\arguments{
  \item{x}{spectrum1 (MSnbase::Spectrum)}

  \item{y}{spectrum2 (MSnbase::Spectrum)}

  \item{fun}{similarity function (must take two spectra and
  ... as arguments)}

  \item{...}{further arguments passed to "fun"}
}
\value{
double, similarity score
}
\description{
calculate similarity between spectra (between their
intensity profile)
}

