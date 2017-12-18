##' Select feature variables to be retained.
##'
##' @title Select feature variables of interest
##'
##' @param object An \code{MSnSet}, \code{MSnExp} or \code{OnDiskMSnExp}.
##' 
##' @param graphics A \code{logical} (default is \code{TRUE})
##'     indicating whether a shiny application should be used if
##'     available. Otherwise, a text menu is used. Ignored if \code{k}
##'     is not missing.
##' 
##' @param fcol A \code{numeric}, \code{logical} or \code{character} of
##'     valid feature variables to be passed directly.
##' 
##' @return For \code{selectFeatureData}: updated object containing only
##'     selected feature variables.
##' 
##' @author Laurent Gatto
##' 
##' @examples
##' 
##' library("pRolocdata")
##' data(hyperLOPIT2015)
##' ## 5 first feature variables
##' x <- selectFeatureData(hyperLOPIT2015, fcol = 1:5)
##' fvarLabels(x)
##' \dontrun{
##' ## select via GUI
##' x <- selectFeatureData(hyperLOPIT2015)
##' fvarLabels(x)
##' }
##'
##' ## Subset the feature data of an OnDiskMSnExp object to the minimal
##' ## required columns
##' f <- system.file("microtofq/MM14.mzML", package = "msdata")
##' od <- readMSData(f, mode = "onDisk")
##'
##' ## what columns do we have?
##' fvarLabels(od)
##'
##' ## Reduce the feature data data.frame to the required columns only
##' od <- selectFeatureData(od, fcol = requiredFvarLabels(class(od)))
##' fvarLabels(od)
selectFeatureData <- function(object,
                              graphics = TRUE,
                              fcol) {
    if (missing(fcol)) {
        if (graphics) {
            if (!requireNamespace("shiny", quietly = TRUE)) {
                warning("The shiny package is required to use the graphical interface.")
                fcol <- .selectTextFeatureData(object)
            } else
                fcol <- .selectShinyFeatureData(object)
        } else fcol <- .selectTextFeatureData(object)
    }
    fData(object) <- fData(object)[, fcol, drop = FALSE]
    if (validObject(object))
        object
}


.selectTextFeatureData <- function(object)
    select.list(fvarLabels(object), multiple=TRUE)


.selectShinyFeatureData <- function(object) {
    sel <- fv <- fvarLabels(object)
    on.exit(return(sel))
    ui <- shiny::fluidPage(
        title = 'Examples of DataTables',
        shiny::sidebarLayout(
                   shiny::sidebarPanel(
                              shiny::actionButton("stop", "Stop app"),
                              shiny::checkboxGroupInput('vars', 'Feature variables',
                                                        as.list(fv), selected = sel)),
                   shiny::mainPanel(shiny::dataTableOutput('fd'))))    
    server <- function(input, output) {
        shiny::observeEvent(input$stop, {
            shiny::stopApp(returnValue = sel)
        })        
        output$fd <- shiny::renderDataTable({
            sel <<- input$vars
            fData(object)[, input$vars, drop = FALSE]
        })
    }
    app <- list(ui=ui, server=server)
    shiny::runApp(app)
}

#' @rdname selectFeatureData
#'
#' @description `requiredFvarLabels` returns a `character` vector with the
#'     required feature data variable names (`fvarLabels`, i.e. the column
#'     names in the `fData` `data.frame`) for the specified object.
#' 
#' @param x `character(1)` specifying the class name for which the required
#'     feature data variable names should be returned.
#'
#' @return For `requiredFvarLabels`: `character` with the required feature
#'     variable names.
#' 
#' @md
requiredFvarLabels <- function(x = c("OnDiskMSnExp", "MSnExp", "MSnSet")) {
    x <- match.arg(x)
    if (x == "OnDiskMSnExp")
        .MSnExpReqFvarLabels
    else character()
}
