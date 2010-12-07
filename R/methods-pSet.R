setMethod("initialize",
          signature(.Object="eSet"),
          function(.Object,
                   assayData,
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData = new("MIAPE"),
                   protocolData = phenoData[,integer(0)],
                   ...) {



          }
