## Count Aggregation by Cell Type ----

#' Within cell type gene-by-subject matrices
#'
#' Given a vector of unique cell types, calculate a gene-by-subject matrix and
#' inter-subject metadata for a differential expression design matrix.
#'
#' @inheritParams aggregateBioVar
#' @export
#'
#' @return List of gene-by-subject and design matrices for each cell type.
#' @seealso \code{\link{summarizedCounts}} for aggregate counts and metadata
#'   \linkS4class{SummarizedExperiment} object.
#'
#' @examples
#' ## Return list of `SummarizedExperiments` with gene-by-subject count matrices
#' ## and subject metadata for each unique `SingleCellExperiment` cell type.
#' countsByCell(
#'     scExp=small_airway,
#'     subjectVar="orig.ident", cellVar="celltype"
#' )
#'
countsByCell <- function(scExp, subjectVar, cellVar) {

    cellTypes <-
        scSubjects(scExp, subjects=subjectVar, cellTypes=cellVar)[["cellTypes"]]
    cellCountList <-
        lapply(cellTypes, function(cellType) {

            ## Subset counts and metadata by cell type into new SCE object.
            cellIndex <- which(scExp[[cellVar]] == cellType)
            cellCounts <- SingleCellExperiment::counts(scExp)[, cellIndex]
            cellMetaData <- SummarizedExperiment::colData(scExp)[cellIndex, ]
            cellAggregate <-
                SingleCellExperiment::SingleCellExperiment(
                    assays=list("counts"=cellCounts),
                    colData=cellMetaData
                )

            ## Aggregate within-cell type gene counts and metadata by subject.
            summarizedCounts(scExp=cellAggregate, subjectVar=subjectVar)
            })

    cellCountList <- stats::setNames(object=cellCountList, nm=cellTypes)
    return(cellCountList)
}

#' Aggregate subject-level biological variation
#'
#' Given an input gene-by-cell count matrix from a
#' \linkS4class{SingleCellExperiment} object, sum within-subject gene counts
#' into an aggregate gene-by-subject count matrix. Column metadata accessed
#' by \link[SummarizedExperiment:colData]{SummarizedExperiment::colData()}
#' are collated by \link{subjectMetaData} to remove variables with inter-cell
#' intrasubject variation, effectively retaining between-subject variation.
#' The summary operations are performed across all cell types and within each
#' cell type. A list of \linkS4class{SummarizedExperiment} objects is returned
#' each with aggregate gene-by-subject count matrix and inter-subject metadata.
#'
#' @param scExp \linkS4class{SingleCellExperiment} object containing
#'   (at minimum) gene counts and column metadata describing sample identifiers
#'   and cell types.
#' @param subjectVar Metadata column name assigning biological sample
#'   identity to aggregate within-subject feature counts.
#' @param cellVar Metadata column name assigning cell type. Used for
#'   aggregating gene-by-subject count matrices by cell type.
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
#'
#' @return List of \linkS4class{SummarizedExperiment} objects with
#'   gene-by-subject count matrices and variable inter-subject column
#'   metadata across and within cell types.
#'
#' @examples
#' ## Aggregate gene-by-subject count matrix and inter-subject metadata
#' aggregateBioVar(
#'     scExp=small_airway,
#'     subjectVar="orig.ident", cellVar="celltype"
#' )
#'
aggregateBioVar <- function(scExp, subjectVar, cellVar) {

    ## Cast single cell input as `SingleCellExperiment` class
    if (!methods::is(scExp, "SingleCellExperiment")) {
        stop(paste0("`scExp` must be an object of SingleCellExperiment class:",
                    "\n\nhttps://bioconductor.org/packages/release/",
                    "bioc/html/SingleCellExperiment.html"))
    }

    ## Notify subsetting for multiple `SingleCellExperiment` assay slots.
    if (length(SummarizedExperiment::assays(scExp)) > 1) {
        assayName <- SummarizedExperiment::assayNames(scExp)[1]
        message(paste("SingleCellExperiment counts from assay:", assayName))
        scExp <- SingleCellExperiment::SingleCellExperiment(
            assay=list("counts"=SummarizedExperiment::assay(
                scExp, assayName
            )),
            colData=SummarizedExperiment::colData(scExp)
        )
    }

    ## Notify metadata variable class to character vector coercion.
    sceClasses <- list(
        "subjectVar"=class(scExp[[subjectVar]]),
        "cellVar"=class(scExp[[cellVar]])
    )
    if (TRUE %in% (!sceClasses %in% "character")) {
        classVars <- names(sceClasses)[!sceClasses %in% "character"]
        varMsg <- vapply(
            X=classVars, FUN.VALUE=character(1),
            FUN=function(var) eval(rlang::sym(var))
        )
        message(
            "Coercing metadata variable to character: ",
            paste(varMsg, collapse=" ")
        )
    }

    ## Combine summarized counts across all cells with aggregation by cell type.
    summarizedAll <- summarizedCounts(scExp=scExp, subjectVar=subjectVar)
    summarizedByCell <-
        countsByCell(scExp=scExp, subjectVar=subjectVar, cellVar=cellVar)
    aggregateSingleCell <- c(list(AllCells=summarizedAll), summarizedByCell)
    return(aggregateSingleCell)
}

