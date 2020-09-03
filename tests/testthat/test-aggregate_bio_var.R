test_that("Cell Type Wrapper", {
    cellSummaries <- countsByCell(
        scExp = small_airway,
        subjectVar = "orig.ident", cellVar = "celltype"
    )
    expect_equivalent(
        names(cellSummaries),
        unique(as.character(small_airway$celltype))
    )
    expect_identical(
        unique(vapply(cellSummaries, class, character(1))),
        "SummarizedExperiment"
    )
})

test_that("Cell Type Wrapper", {
    library(SummarizedExperiment)
    library(SingleCellExperiment)

    ## Error for non-`SingleCellExperiment` object input for `scExp`.
    seWarning <-
        SummarizedExperiment(
            assays = list("counts" = assay(small_airway, "counts")),
            colData = colData(small_airway)
        )
    expect_error(
        aggregateBioVar(
            scExp = seWarning, subjectVar = "orig.ident", cellVar = "celltype"
        )
    )

    ## Message to subset first assay slot if multiple slots are present.
    scExpMultiAssay <- SingleCellExperiment(
        assay = list(
            "RNA" = assay(small_airway, "counts"),
            "counts" = assay(small_airway, "counts")
        ),
        colData = colData(small_airway)
    )
    expect_message(
        aggregateCounts <- aggregateBioVar(
            scExp = scExpMultiAssay,
            subjectVar = "orig.ident", cellVar = "celltype"
        ), regexp = "counts from assay: RNA"
    )

    ## Message for character vector coercion of metadata variables.
    expect_message(
        aggregateBioVar(
            scExp = small_airway,
            subjectVar = "orig.ident", cellVar = "celltype"
        ), regexp = "Coercing metadata variable to character: celltype"
    )

    ## Expected return list output with cell type names.
    expect_is(aggregateCounts, "list")
    expect_equivalent(
        c("AllCells", unique(as.character(small_airway$celltype))),
        names(aggregateCounts)
    )
})

