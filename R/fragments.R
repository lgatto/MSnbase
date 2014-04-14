calculateFragments <- function(sequence, modifiedCystein=TRUE) {
  ## constants; source wikipedia
  proton <- 1.007276466
  water <- 1.0078250321 * 2 + 15.9949146221

  aa <- .get.amino.acids()
  aamass <- setNames(aa$ResidueMass, aa$AA)

  if (modifiedCystein) {
    aamass["C"] <- 160.03065
  }

  fragment.seq <- strsplit(sequence, "")[[1]]
  n <- length(fragment.seq)

  b <- setNames(cumsum(aamass[fragment.seq]) + proton,
                paste0("b", 1:n))
  y <- setNames(cumsum(aamass[rev(fragment.seq)]) + water + proton,
                paste0("y", 1:n))

  fragment.seq <- c(substring(sequence, rep(1, n), 1:n),
                    rev(substring(sequence, 1:n, rep(n, n))))
  fragment.str <- c(names(b), names(y))
  mass <- c(b, y)

  o <- order(mass)

  return(data.frame(mass=mass[o],
                    fragment.str=fragment.str[o],
                    fragment.seq=fragment.seq[o],
                    stringsAsFactors=FALSE))
}

