#' SingleCellExperiment subject values
#'
#' Extract unique values from \linkS4class{SingleCellExperiment}
#' column (i.e. cell) metadata. Used to determine subject and cell type values.
#'
#' @param ... Named metadata variables for subjects and cell types.
#' @inheritParams aggregateBioVar
#' @export
#'
#' @return List of character vectors with unique values from
#' \linkS4class{SingleCellExperiment} column metadata variables.
#'
#' @examples
#' ## Examples of metadata column variable names.
#' names(SummarizedExperiment::colData(small_airway))
#'
#' ## Return list of subject and cell type values from experiment metadata.
#' scSubjects(scExp=small_airway, subjects="orig.ident", cellTypes="celltype")
#'
scSubjects <- function(scExp, ...) {
    metavars <- list(...)
    stopifnot(!is.null(names(metavars)))
    metavar_names <- vapply(
        X=names(metavars), FUN.VALUE=logical(1),
        FUN=function(metavar_name) nchar(metavar_name) > 0
    )
    if (FALSE %in% metavar_names) {
        stop("Names required for: ", metavars[which(metavar_names == FALSE)])
    }
    meta_values <- lapply(metavars, function(var) {
        unique(as.character(scExp[[var]]))
    })
    names(meta_values) <- names(metavars)
    return(meta_values)
}

#' Gene-by-subject count matrix
#'
#' Convert gene-by-cell count matrix to gene-by-subject count matrix.
#' Row sums are calculated for each feature (i.e. gene) across cells by subject.
#'
#' @inheritParams aggregateBioVar
#' @importFrom rlang !! :=
#' @export
#'
#' @return S4 DataFrame of gene-by-subject count sums.
#' @seealso \code{\link{scSubjects}} for `subjects` values.
#'
#' @examples
#' ## Return cell count matrix aggregated by subject.
#' countsBySubject(scExp=small_airway, subjectVar="orig.ident")
#'
countsBySubject <- function(scExp, subjectVar) {

    subjects <- scSubjects(scExp=scExp, subjects=subjectVar)[["subjects"]]

    subjectCounts <-
        lapply(subjects, function(subjectID) {
            subjectMatches <- which(scExp[[subjectVar]] == subjectID)
            if (length(subjectMatches) == 1) {
                subjectSums <-
                    SingleCellExperiment::counts(scExp)[, subjectMatches]
            } else {
                sce <- SingleCellExperiment::counts(scExp)[, subjectMatches]
                subjectSums <- Matrix::rowSums(x=sce)
            }
            tibble::tibble(!!subjectID := unname(subjectSums))
        })

    subjectCounts <- S4Vectors::DataFrame(subjectCounts)
    rownames(subjectCounts) <- rownames(scExp)
    return(subjectCounts)
}

#' Collate metadata variation
#'
#' Identify single cell experiment metadata variables that are identical
#' within subject (e.g. genotype, treatment, cell line, sample preparation).
#' Effectively excludes metadata variables containing between cell variation.
#' Used as the design matrix for differential expression analysis.
#'
#' @inheritParams aggregateBioVar
#' @export
#'
#' @return Tibble data frame of metadata variables without intrasubject
#'   variation from single cell experiment metadata. Rows correspond to
#'   aggregated cells (i.e. subject / biological replicate) and columns to
#'   metadata attribute variables (e.g. genotype, treatment, cell line).
#'
#' @examples
#' ## Return experiment metadata sans intrasubject variation.
#' subjectMetaData(scExp=small_airway, subjectVar="orig.ident")
#'
subjectMetaData <- function(scExp, subjectVar) {

    subjects <- scSubjects(scExp=scExp, subjects=subjectVar)[["subjects"]]
    indexCellBySubject <- match(subjects, scExp[[subjectVar]])

    ## For each metadata variable, calculate a contingency table of variable
    ## levels by subject. Retain variables with a single intrasubject level.
    indexMetaVar <-
        vapply(
            X=SummarizedExperiment::colData(scExp), FUN.VALUE=logical(1),
            FUN=function(metaVar) {
                crosstab <- table(metaVar, scExp[[subjectVar]])
                ifelse(
                    test=TRUE %in% (colSums(crosstab > 0) > 1),
                    yes=FALSE, no=TRUE)
        })

    ## Subset and remove barcode rowname carryover from metadata rownames.
    subjectVariation <-
        SummarizedExperiment::colData(scExp)[indexCellBySubject, indexMetaVar]
    rownames(subjectVariation) <- subjectVariation[[subjectVar]]
    return(subjectVariation)
}

#' Aggregate feature counts and metadata by subject
#'
#' Given an input sparse count matrix and corresponding column metadata,
#' aggregate gene counts by subject level. Metadata variables with only
#' inter-subject variation are retained; any variables with cell-level variation
#' within a subject are dropped (e.g. feature / RNA count by cell).
#'
#' @inheritParams aggregateBioVar
#' @export
#'
#' @return \linkS4class{SummarizedExperiment} object with feature counts
#'   aggregated by subject and summarized inter-subject metadata.
#'
#' @examples
#' ## Construct SummarizedExperiment object with gene-by-subject count matrix
#' ## and column metadata summarized to exclude intrasubject variation.
#' ## See `SummarizedExperiment` accessor functions `assay()` and `colData()`
#' ## to access the count matrix and column metadata for downstream analyses.
#' summarizedCounts(scExp=small_airway, subjectVar="orig.ident")
#'
summarizedCounts <- function(scExp, subjectVar) {
    subjectCounts <- countsBySubject(scExp=scExp, subjectVar=subjectVar)
    summarizedMetadata <- subjectMetaData(scExp=scExp, subjectVar=subjectVar)
    countSummary <-
        SummarizedExperiment::SummarizedExperiment(
            assays=list(counts=subjectCounts),
            colData=summarizedMetadata
        )
    return(countSummary)
}

