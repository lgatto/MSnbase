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
        
        # Extract metadata
        metaData <- list(
            processingData = processingData(from),
            experimentData = experimentData(from),
            protocolData = protocolData(from),
            annotation = annotation(from),
            qual = qual(from))

        # Generate SE
        SummarizedExperiment::SummarizedExperiment(
            assays = as.matrix(raw),
            rowData = rowData,
            colData = colData,
            metadata = metaData)
})

setAs("SummarizedExperiment",
    "MSnSet",
    function(from) {
        if(!requireNamespace("SummarizedExperiment")) {
            stop("The SummarizedExperiment package is required",
                "for SummarizedExperiment to MSnSet conversion.")
        }
        
        # Extract expression, feature and pheno data
        raw <- SummarizedExperiment::assay(from)
        featData <- data.frame(
            SummarizedExperiment::rowData(from), 
            row.names = names(from))
        phenoData <- data.frame(SummarizedExperiment::colData(from))
        
        # Extract metadata
        processingData = metadata(from)$processingData
        experimentData = metadata(from)$experimentData
        protocolData = metadata(from)$protocolData
        annotation = metadata(from)$annotation
        qual = metadata(from)$qual

        # Generate MSnSet
        msnset <- MSnSet(exprs = as.matrix(raw),
            pData = AnnotatedDataFrame(phenoData),
            fData = AnnotatedDataFrame(featData))
        if(!is.null(experimentData)) 
            msnset@experimentData <- experimentData
        if(!is.null(processingData)) 
            msnset@processingData <- processingData
        if(!is.null(protocolData)) 
            msnset@protocolData <- protocolData
        if(!is.null(annotation)) 
            msnset@annotation <- annotation
        if(!is.null(qual)) 
            msnset@qual <- qual
        return(msnset)
})
