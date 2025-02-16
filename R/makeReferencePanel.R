#' Function for assembling a SummarizedExperiment object of reference panel in \strong{hg19} for
#' cell type deconvolution (which is used in \code{decemedip} function)
#'
#' @param row_ranges A \code{GRanges} object that contains the genomic coordinates
#' of \emph{reference regions/sites}.
#' @param X A matrix that contains the beta values of reference regions. \strong{Each row
#' is a region and each column is a cell type.} If \code{X} has row names or column names,
#' the output \code{SummarizedExperiment} object will inherit them.
#' @param cpg_coords A \code{GRanges} object that contains genomic coordinates of
#' all CpGs in the genome. A ready-to-use CpG list for hg19 is available to download at
#' \url{https://github.com/nshen7/decemedip-scripts/blob/main/cpg_coords/hg19.cpg.coords.rda}.
#' It is used for generating coloum `n_cpg_100bp` in the reference panel, which
#' represents CpG density around the reference CpG.
#' @param col_names An atomic vector of strings that indicates the column names, i.e.,
#' names of the cell types. Default is NULL. If not NULL, the column names of \code{X}
#' will be overwritten by this argument.
#' @param row_names An atomic vector of strings that indicates the row names, i.e.,
#' names of the reference regions. Default is NULL. If not NULL, the row names of \code{X}
#' will be overwritten by this argument.
#' @param col_data A \code{DataFrame} object that contains metadata for columns (i.e.,
#' cell types) if specified. Each row in \code{col_data} should contain info of a cell type
#' in \code{X}. If input is a non-\code{DataFrame} object, it will be converted
#' to a \code{DataFrame}. Default is NULL.
#' @param row_data A \code{DataFrame} object that contains metadata for row (i.e.,
#' reference regions) if specified. Each row in \code{row_data} should contain info of a
#' reference region in \code{X}. If input is a non-\code{DataFrame} object, it will
#' be converted to a \code{DataFrame}. Default is NULL.
#'
#' @importFrom methods is
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom GenomicRanges granges resize
#' @importFrom GenomicRanges countOverlaps
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom GenomicRanges GRanges
#'
#' @return An \code{SummarizedExperiment} object with each row represents a reference region
#' and an assay named 'X' that stores the beta values of reference regions.
#' @export
#'
#' @examples
#' gr <- GenomicRanges::GRanges(
#'   seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr3")),      # Chromosome names
#'   ranges = IRanges::IRanges(start = c(101123, 203456, 987654), # Start positions
#'                             end = c(101125, 203458, 987656)),  # End positions
#'   cpg_id = c("CpG_001", "CpG_002", "CpG_003")      # CpG site IDs
#' )
#' makeReferencePanel(
#'   row_ranges = gr,
#'   X = matrix(runif(6), nrow = 3)
#' )
#'
makeReferencePanel <- function(
    row_ranges,
    X,
    cpg_coords,
    col_names = NULL,
    row_names = NULL,
    col_data = NULL,
    row_data = NULL
) {

  if (length(row_ranges) != nrow(X)) stop('Length of `row_ranges` is not equal to number of rows in `X`!')

  if (!is.null(col_data) & !methods::is(col_data, 'DataFrame')) {
    message('col_data is being converted to a DataFrame object.')
    col_data <- DataFrame(col_data)
  }
  if (!is.null(row_data) & !methods::is(row_data, 'DataFrame')) {
    message('row_data is being converted to a DataFrame object.')
    row_data <- DataFrame(row_data)
  }
  if (!is.null(col_names))
    stopifnot('The length of `col_names` not equal to number of columns in X.' = length(col_names) == ncol(X))
  if (!is.null(row_names))
    stopifnot('The length of `row_names` not equal to number of rows in X.' = length(row_names) == nrow(X))
  if (!is.null(col_data))
    stopifnot('The number of rows in col_data not equal to number of columns in X.' = nrow(col_data) == ncol(X))
  if (!is.null(row_data))
    stopifnot('The number of rows in row_data not equal to number of rows in X.' = nrow(row_data) == nrow(X))

  if (!is.null(col_names)) colnames(X) <- col_names
  if (!is.null(row_names)) rownames(X) <- row_names

  se <- SummarizedExperiment::SummarizedExperiment(rowRanges = row_ranges, assays = list(beta = X))

  if (!is.null(col_data)) {
    if (!is.null(col_names)) rownames(col_data) <- col_names
    colData(se) <- cbind(colData(se), col_data)
  }
  if (!is.null(row_data)) {
    if (!is.null(row_names)) rownames(row_data) <- row_names
    rowData(se) <- cbind(rowData(se), row_data)
  }

  ## Include CpG density as covariate
  rowData(se)$n_cpgs_100bp <- IRanges::countOverlaps(
    GenomicRanges::granges(se) |>
      GenomicRanges::resize(width = 100, fix = 'center'),
    cpg_coords)

  return(se)
}
