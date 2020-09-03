#' Small Airway `SingleCellExperiment` object
#'
#' Small airway epithelium single-cell RNA sequencing data subset combined from
#' 7 porcine individuals. Genotypes represent wildtype (Genotype `WT`; n=3)
#' and CFTR-knockout subjects (Genotype `CFTRKO`; n=4) expressing a
#' cystic fibrosis phenotype.\cr
#' \cr
#' To access cell counts and column metadata:
#' \itemize{
#'   \item counts: `SummarizedExperiment::assay(small_airway, "counts")`
#'   \item metadata: `SummarizedExperiment::colData(small_airway)`
#'   }
#' Cell types include:
#' \itemize{
#'   \item Secretory cell
#'   \item Endothelial cell
#'   \item Immune cell
#' }
#'
#' @format A \linkS4class{SingleCellExperiment} object with feature counts of
#' 1311 genes (rows) from 2687 individual cells (columns).
#' Features were subset by gene ontology annotation "ion transport" with
#' included child terms ("anion transport", "cation transport", etc.).
#' Includes a count matrix accessed by \link[SummarizedExperiment]{assay} and
#' column metadata accessed by \link[SummarizedExperiment]{colData}.
#' Metadata variables are defined as:
#' \describe{
#'   \item{orig.ident}{Biological sample identifier (i.e. subject)}
#'   \item{nCount_RNA}{Number of UMIs (i.e. depth/size factor)}
#'   \item{nFeature_RNA}{Features with non-zero counts
#'     (i.e. identified features)}
#'   \item{Genotype}{Sample Genotype ('CF' == cystic fibrosis)
#'     \itemize{
#'       \item WT = non-CF pig
#'       \item CFTRKO = CFTR knockout}
#'     }
#'   \item{Sample}{Airway region of sample}
#'   \item{celltype}{Cell type labels}
#' }
#'
"small_airway"

