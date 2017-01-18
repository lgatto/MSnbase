## These data should eventually be extracted in mzR and added to the
## default header

injectionTimeFromFile1 <- function(f) {
    xvals <- vector("list", length = length(f))
    doc <- XML::xmlInternalTreeParse(f)
    namespaceDef <- XML::getDefaultNamespace(doc)
    if (length(namespaceDef) == 0)
        return(NA)
    ns <- c(x = namespaceDef[[1]]$uri)
    x <- XML::xpathApply(doc, "//x:cvParam[@accession = 'MS:1000927']",
                         namespaces = ns, xmlAttrs)
    xvals <- sapply(x, "[", "value")
    as.numeric(unlist(xvals))
}
