## suggested use: keep only non-empty slots - see .reduce_list in
## utils.R

setAs("AnnotatedDataFrame", "list",
      function (from) as.list(from@data))

setAs("MIAxE", "list", 
      function (from) {
          nms <- slotNames(from)
          nms <- setdiff(nms, ".__classVersion__")
          ans <- vector("list", length = length(nms))
          names(ans) <- nms
          for (k in nms)
              ans[[k]] <- slot(from, k)
          ans
      })

setAs("MSnProcess", "list", 
      function (from) {
          nms <- slotNames(from)
          nms <- setdiff(nms, ".__classVersion__")
          ans <- vector("list", length = length(nms))
          names(ans) <- nms
          for (k in nms)
              ans[[k]] <- slot(from, k)
          ans
      })

