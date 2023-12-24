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

lockEnvironment(.MSnbaseEnv,bindings=TRUE)
