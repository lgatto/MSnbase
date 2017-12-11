setAs("MSnSet",
    "SummarizedExperiment",
    function(from) {
        if(!requireNamespace("SummarizedExperiment")) {
            stop("The SummarizedExperiment package is required",
                "for MSnSet to SummarizedExperiment conversion.")
        }

        # Extract assay, rowData and colData
        raw <- exprs(from)
        rowData <- fData(from)
        colData <- pData(from)
        
        # Generate SE
        SummarizedExperiment(assays = as.matrix(raw),
            rowData = rowData,
            colData = colData)
})

setAs("SummarizedExperiment",
    "MSnSet",
    function(from) {
        if(!requireNamespace("SummarizedExperiment")) {
            stop("The SummarizedExperiment package is required",
                "for SummarizedExperiment to MSnSet conversion.")
        }
        
        # Extract expression, feature and pheno data
        raw <- assay(from)
        featData <- data.frame(rowData(from), row.names = names(from))
        phenoData <- data.frame(colData(from))

        # Generate MSnSet
        MSnSet(exprs = as.matrix(raw),
            pData = AnnotatedDataFrame(phenoData),
            fData = AnnotatedDataFrame(featData))
})