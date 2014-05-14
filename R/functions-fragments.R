#' calculate fragments from a peptide sequence
#' @param sequence character vector of length 1
#' @param type could be c("a", "b", "c", "x", "y", "z")
#' @param z charge
#' @param modifications a named (amino acid one-letter-code; upper case) vector
#' of modified mass (default: Carbamidomethyl (C) replaces Cystein: 160.030649).
#' @param verbose verbose output?
#' @noRd
.calculateFragments <- function(sequence, type=c("b", "y"), z=1,
                                modifications=c(C=160.030649), verbose=TRUE) {
  type <- match.arg(type, choices=c("a", "b", "c", "x", "y", "z"), several.ok=TRUE)
  type <- sort(type)
  ## constants
  mass <- get.atomic.mass()
  ## according to Table 1 of:
  ## Johnson, R. S., Martin, S. A., Biemann, K., Stults, J. T., and
  ## Watson, J. T. (1987).
  ## Novel fragmentation process of peptides by collision-induced
  ## decomposition in a tandem mass spectrometer: differentiation of leucine
  ## and isoleucine.
  ## Analytical Chemistry, 59(21), 2621-2625.
  ## http://dx.doi.org/10.1021/ac00148a019
  ##
  ## a proton (H+) is added later
  ## (after calculation of the different charge states)
  add <- c(a=-(mass["C"]+mass["O"]),            # + H - CO
           b=0,                                 # + H
           c=mass["N"]+3*mass["H"],             # + H + NH3
           x=mass["C"]+2*mass["O"],             # + CO + OH
           y=2*mass["H"]+mass["O"],             # + H2 + OH
           z=-(mass["N"]+mass["H"])+mass["O"])  # + NH + OH

  aa <- .get.amino.acids()
  aamass <- setNames(aa$ResidueMass, aa$AA)

  ## replace default mass by modifications
  if (length(modifications)) {
    aamass[names(modifications)] <- modifications

  }

  if (verbose) {
    if(length(modifications)) {
      mods <- paste0(names(modifications), "=", modifications, collapse=", ")
    } else {
      mods <- "None"
    }
    message("Modifications used: ", mods)
  }

  ## split peptide sequence into aa
  fragment.seq <- strsplit(sequence, "")[[1]]

  ## calculate cumulative mass starting at the amino-terminus (for a, b, c)
  amz <- cumsum(aamass[fragment.seq])
  ## calculate cumulative mass starting at the carboxyl-terminus (for x, y, z)
  cmz <- cumsum(aamass[rev(fragment.seq)])

  ## calculate fragment mass (amino-terminus)
  tn <- length(amz)
  atype <- c("a", "b", "c") %in% type
  nat <- sum(atype)
  amz <- rep(amz, nat) + rep(add[1:3][atype], each=tn)

  ### calculate fragment mass (carboxyl-terminus)
  ctype <- c("x", "y", "z") %in% type
  nct <- sum(ctype)
  cmz <- rep(cmz, nct) + rep(add[4:6][ctype], each=tn)

  ## devide by charge
  zn <- length(z)
  amz <- rep(amz, each=zn)/z
  cmz <- rep(cmz, each=zn)/z

  ## add protons (H+)
  amz <- amz + mass["p"]
  cmz <- cmz + mass["p"]

  ## fragment seq (amino-terminus)
  fn <- length(fragment.seq)
  aseq <- rep(rep(substring(sequence, rep(1, fn), 1:fn), each=zn), nat)

  ## fragment seq (carboxyl-terminus)
  cseq <- rep(rep(rev(substring(sequence, 1:fn, rep(fn, fn))), each=zn), nct)

  ## fragment str (amino-terminus)
  atype <- rep(c("a", "b", "c")[atype], each=tn*zn)
  pos <- rep(1:tn, each=zn)
  aion <- paste0(atype, pos)

  ## fragment str (carboxyl-terminus)
  ctype <- rep(c("x", "y", "z")[ctype], each=tn*zn)
  cion <- paste0(ctype, pos)

  df <- data.frame(mz=c(amz, cmz),
                   ion=c(aion, cion),
                   type=c(atype, ctype),
                   pos=pos,
                   z=z,
                   seq=c(aseq, cseq),
                   stringsAsFactors=FALSE)
  return(df)
}

