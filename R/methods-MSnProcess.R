##################################################################
## Methods for MSnProcess class
setMethod("show","MSnProcess",function(object) show.MSnProcess(object))

setMethod("fileNames",
          signature(object="MSnProcess"),
          function(object) object@files)
