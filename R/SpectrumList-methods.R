#' @include SpectrumList.R

setMethod("show", "SpectrumList", function(object) {
    .show_SpectrumList(object, margin = "  ", print.classinfo = TRUE)
})

