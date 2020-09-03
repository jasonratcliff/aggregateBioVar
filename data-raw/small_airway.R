# Load small airway subset Seurat object
# - Includes 1311 genes from 2687 cells (secretory, endothelial, immune)
# - 3 CF and 4 non-CF genotype biological (porcine) replicates
load("data-raw/ASE_integrated_sub.Rda")

# Create SingleCellExperiment object from RNA assay counts and Seurat metadata.
# - Count dgCMatrix from Seurat object RNA Assay
# - Metadata data frame from Seurat object
small_airway <-
    SingleCellExperiment::SingleCellExperiment(
        assays = list("counts" = Seurat::GetAssay(object = ASE.integrated.sub,
                                                  "RNA")@counts),
        colData = ASE.integrated.sub@meta.data
    )

# Replace "Sample" column name with "Region"
colDataNames <- names(SummarizedExperiment::colData(small_airway))
names(SummarizedExperiment::colData(small_airway)) <-
    gsub(pattern = "^Sample$", replacement = "Region", x = colDataNames)

usethis::use_data(small_airway)
