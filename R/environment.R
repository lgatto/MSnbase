.MSnbaseEnv <- new.env(parent = emptyenv(), hash = TRUE)

## As discussed in issue #163 for details, the random errors we see
## (see issue #138) seem to come (partially at least) from using new
## in the prototype. As a result, these will be setn (and tested in
## validity methods) outside of the prototype. The vector below stores
## the respective class versions. When a class doesn't have one, the
## version should be defined as NA_character_.

ClassVersions <- c(
    Spectrum = "0.4.0",
    Spectrum1 = "0.2.0",
    Spectrum2 = "0.2.0")

assign("ClassVersions", ClassVersions, envir = .MSnbaseEnv)

assign("amino.acids",
       data.frame(AA = c("peg","A","R","N","D","C","E",
                    "Q","G","H","I","L", "K","M","F",
                    "P","S","T","W","Y","V"),
                  ResidueMass = c(44.00000,
                    71.03711,  156.10111, 114.04293, 115.02694,
                    103.00919, 129.04259, 128.05858, 57.02146,
                    137.05891, 113.08406, 113.08406, 128.09496,
                    131.04049, 147.06841, 97.05276,  87.03203,
                    101.04768, 186.07931, 163.06333, 99.06841),
                  Abbrev3 = c(NA,
                    "Ala", "Arg", "Asn", "Asp", "Cys",
                    "Glu", "Gln", "Gly", "His", "Ile",
                    "Leu", "Lys", "Met", "Phe", "Pro",
                    "Ser", "Thr", "Trp", "Tyr", "Val"),
                  ImmoniumIonMass = c(NA,
                    44.05003,  129.11400, 87.05584,  88.03986,  76.02210,
                    102.05550, 101.07150, 30.03438,  110.07180, 86.09698,
                    86.09698,  101.10790, 104.05340, 120.08130, 70.06568,
                    60.04494,  74.06059,  159.09220, 136.07620, 72.08133),
                  Name = c("Polyethylene glycol",
                    "Alanine",    "Arginine",      "Asparagine", "Aspartic acid",
                    "Cysteine",   "Glutamic acid", "Glutamine",  "Glycine",
                    "Histidine",  "Isoleucine",    "Leucine",    "Lysine",
                    "Methionine", "Phenylalanine", "Proline",    "Serine",
                    "Threonine",  "Tryptophan",    "Tyrosine",   "Valine"),
                  ## The hydrophobicity values are from JACS, 1962, 84: 4240-4246. (C. Tanford)
                  Hydrophobicity = c(NA, 0.62, -2.53, -0.78, -0.9, 0.29,
                    -0.74, -0.85, 0.48, -0.4, 1.38, 1.06, -1.5, 0.64, 1.19,
                    0.12, -0.18, -0.05, 0.81, 0.26, 1.08),
                  ## The hydrophilicity values are from PNAS, 1981, 78:3824-3828 (T.P.Hopp & K.R.Woods)
                  Hydrophilicity = c(NA, -0.5, 3, 0.2, 3, -1, 3, 0.2, 0,
                    -0.5, -1.8, -1.8, 3, -1.3, -2.5, 0, 0.3, -0.4, -3.4,
                    -2.3, -1.5),
                  SideChainMass = c(NA, 15, 101, 58, 59, 47, 73, 72,
                    1, 82, 57, 57, 73, 75, 91, 42, 31, 45, 130, 107, 43),
                  ## CRC Handbook of Chemistry and Physics, 66th ed., CRC Press, Boca Raton, Florida (1985).
                  ## R.M.C. Dawson, D.C. Elliott, W.H. Elliott, K.M. Jones, Data for Biochemical Research 3rd ed., Clarendon Press Oxford
                  pK1 = c(NA, 2.35, 2.18, 2.18, 1.88, 1.71, 2.19,
                    2.17, 2.34, 1.78, 2.32, 2.36, 2.2, 2.28, 2.58,
                    1.99, 2.21, 2.15, 2.38, 2.2, 2.29),
                  pK2 = c(NA, 9.87, 9.09, 9.09, 9.6, 10.78, 9.67, 9.13,
                    9.6, 8.97, 9.76, 9.6, 8.9, 9.21, 9.24, 10.6, 9.15,
                    9.12, 9.39, 9.11, 9.74),
                  pI = c(NA, 6.11, 10.76, 10.76, 2.98, 5.02, 3.08,
                    5.65, 6.06, 7.64, 6.04, 6.04, 9.47, 5.74, 5.91,
                    6.3, 5.68, 5.6, 5.88, 5.63, 6.02)),
       envir = .MSnbaseEnv)

assign("atomic.mass",
          ## taken from:
          ## http://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
          c(H=1.007825,
            C=12,
            N=14.003074,
            O=15.994915,
            p=1.007276),
       envir = .MSnbaseEnv)

lockEnvironment(.MSnbaseEnv,bindings=TRUE)
