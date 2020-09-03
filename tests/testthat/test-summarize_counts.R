context("SingleCellExperiment metadata")

test_that("Unique subject and cell type values", {
    ## `orig.ident` == character vector, `celltype` == factor vector
    assayInfo <-
        scSubjects(
            scExp = small_airway,
            subjects = "orig.ident", cellTypes = "celltype"
        )
    expect_equal(
        assayInfo$subjects, unique(small_airway$orig.ident)
    )
    expect_equal(
        assayInfo$cellTypes, as.character(unique(small_airway$celltype))
    )
    ## stopifnot() error for null names
    expect_error(
        scSubjects(scExp = small_airway, "orig.ident", "celltype")
    )
    ## stop() error for missing name in dot param input
    expect_error(
        scSubjects(scExp = small_airway, "orig.ident", cellTypes = "celltype")
    )

})

test_that("Feature count aggregation by subject", {
    subjectCounts <-
        countsBySubject(scExp = small_airway, subjectVar = "orig.ident")
    expect_s4_class(subjectCounts, "DataFrame")
    expect_identical(rownames(subjectCounts), rownames(small_airway))
    expect_error(
        countsBySubject(scExp = small_airway, "orig.ident", "celltype")
    )

    ## From immune cell subset, assign list of cell indexes by subject.
    ic_cells <- small_airway[, grep(pattern = "Immune cell",
                                    x = as.character(small_airway$celltype))]
    subject_index <-
        lapply(
            X = unique(small_airway$orig.ident),
            FUN = function(subject) {
                grep(pattern = subject, x = ic_cells$orig.ident)
        })

    ## Test for single subject matches to prevent rowSums() dimension error.
    subject_sample <- sample(length(subject_index), size = 2)
    subject_subset <-
        subset(subject_index, ifelse(seq_along(subject_index) %in%
                                        subject_sample, FALSE, TRUE))
    subset_index <- c(
        unlist(subject_subset),
        subject_index[subject_sample][[1]][1],
        subject_index[subject_sample][[2]][1]
    )
    cell_subset <- ic_cells[, subset_index]
    expect_gte(
        length(which(table(cell_subset$orig.ident) %in% 1 == TRUE)), 2
    )
    expect_silent(countsBySubject(cell_subset, subjectVar = "orig.ident"))
})

test_that("Subject metadata summary", {
    summarizedMetadata <-
        subjectMetaData(scExp = small_airway, subjectVar = "orig.ident")
    expect_s4_class(summarizedMetadata, "DataFrame")
    expect_equivalent(
        summarizedMetadata$orig.ident,
        unique(small_airway$orig.ident)
    )
    expect_equivalent(
        unique(summarizedMetadata$Genotype),
        unique(small_airway$Genotype)
    )
    expect_equivalent(
        rownames(summarizedMetadata),
        unique(small_airway[["orig.ident"]])
    )
})

test_that("Aggregate SummarizedExperiment output", {
     summarizedSubjects <-
        summarizedCounts(scExp = small_airway, subjectVar = "orig.ident")
    expect_identical(
        SummarizedExperiment::assay(summarizedSubjects, "counts"),
        countsBySubject(scExp = small_airway, subjectVar = "orig.ident")
    )
    expect_identical(
        SummarizedExperiment::colData(summarizedSubjects),
        subjectMetaData(scExp = small_airway, subjectVar = "orig.ident")
    )
    expect_s4_class(summarizedSubjects, "SummarizedExperiment")
})

