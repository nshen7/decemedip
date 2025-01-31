#' Obtain read counts of regions of interest from MeDIP-seq bam files
#'
#' @param sample_bam_files An atomic vector that contains directory paths of the MeDIP-seq
#' sorted bam files.
#' @param sample_names An atomic vector of strings that indicates the sample names.
#' Default is NULL. If not NULL, please make sure that \code{sample_names} corresponds to
#' elements of \code{sample_bam_files}.
#' @param sample_paired A logic value that indicates sample_paired-end reads (if TRUE) or single-end
#' reads (if FALSE).
#' @param roi A \code{GRanges} object that contains the genomic coordinates of the region
#' of interest (ROI). Default is a \code{GRanges} object that contains the default reference set of
#' cell type-specific CpGs.
#' @param col_data A \code{DataFrame} object that contains metadata for columns (i.e.,
#' samples) if specified. Default is NULL. If not NULL, please make sure that
#' rows of \code{col_data} corresponds to elements of sample_bam_files. If input is a
#' non-\code{DataFrame} object, it will be converted to a \code{DataFrame}.
#' @param row_data A \code{DataFrame} object that contains metadata for rows (i.e.,
#' genomic regions) if specified. Default is NULL. If not NULL, please make sure that
#' rows of \code{row_data} corresponds to elements of sample_bam_files. If input is a
#' non-\code{DataFrame} object, it will be converted to a \code{DataFrame}.
#' @param bs_genome A character value that indicates the reference genome name as
#' defined by \code{BSgenome} package. Default is 'BSgenome.Hsapiens.UCSC.hg19'.
#' @param ... Additional arguments passed into \code{MEDIPS::MEDIPS.createROIset}
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom methods is
#' @importFrom S4Vectors DataFrame
#' @importFrom MEDIPS MEDIPS.createROIset
#' @importFrom MEDIPS genome_count
#' @importFrom purrr map_dfc
#' @importFrom SummarizedExperiment colData<-
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom GenomicRanges GRanges
#'
#' @return An object of class \code{SummarizedExperiment} with read count matrix stored
#' as an assay named 'counts' (can be extracted using \code{SummarizedExperiment::assays})
#' @export
#'
#' @examples
#' ?getRoiReadCount
#'
getRoiReadCount <- function(
    sample_bam_files,
    sample_names,
    sample_paired,
    roi = granges(hg19.ref.cts.se),
    col_data = NULL,
    row_data = NULL,
    bs_genome = 'BSgenome.Hsapiens.UCSC.hg19',
    ...
  ) {


  if (!methods::is(roi, 'GRanges'))
    stop('The roi is not a GRanges object.')
  if (!is.null(col_data) & !methods::is(col_data, 'DataFrame')) {
    message('col_data is being converted to a DataFrame object.')
    col_data <- DataFrame(col_data)
  }
  if (!is.null(row_data) & !methods::is(row_data, 'DataFrame')) {
    message('row_data is being converted to a DataFrame object.')
    row_data <- DataFrame(row_data)
  }
  if (!is.null(sample_names))
    stopifnot('The length of sample_names not equal to length of sample_bam_files.' = length(sample_names) == length(sample_bam_files))
  if (!is.null(col_data))
    stopifnot('The number of rows in col_data not equal to length of sample_bam_files.' = nrow(col_data) == length(sample_bam_files))
  if (!is.null(row_data))
    stopifnot('The number of rows in row_data not equal to number of regions.' = nrow(row_data) == length(roi))

  .getRoiReadCountPerSample <- function(bam_dir, ...) {

    roi.df <- as.data.frame(roi, row.names = NULL)[, 1:3] %>%
      dplyr::rename('chr' = 'seqnames') |>
      dplyr::mutate(name = 1:length(roi))
    count_per <- MEDIPS::MEDIPS.createROIset(
      file = bam_dir,
      ROI = roi.df,
      BSgenome = bs_genome,
      sample_name = sample_names,
      paired = sample_paired,
      uniq = 0,
      ...
    ) |> MEDIPS::genome_count()

    return(count_per)
  }

  ## TODO: remove message of column name conversion
  counts.df <- purrr::map_dfc(sample_bam_files, .getRoiReadCountPerSample, .progress = TRUE) |> as.data.frame() |> as.matrix()
  if (!is.null(sample_names)) colnames(counts.df) <- sample_names
  # if (!is.null(row_data)) if (!is.null(rownames(row_data))) rownames(counts.df) <- rownames(row_data)

  se <- SummarizedExperiment::SummarizedExperiment(rowRanges = roi, assays = list(counts = counts.df))

  if (!is.null(col_data)) {
    if (!is.null(sample_names)) rownames(col_data) <- sample_names
    colData(se) <- col_data
  }
  if (!is.null(row_data)) {
    rowData(se) <- row_data
  }

  return(se)
}
