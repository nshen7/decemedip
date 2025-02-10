#' Hg19 genomic information of anchor CpGs (i.e., all-tissue unmethylated/methylation
#' probes) inferred from DNA methylation atlas published by Moss 2018 Nat. Commun.
#' (https://www.nature.com/articles/s41467-018-07466-6). Used as default in \code{\link{decemedip}}.
#'
#'
#' @description This dataset represents a GRanges object that contains the collection
#' of Illumina HumanMethylation450K probes that have methylation level less than 0.1 or
#' greater than 0.9 in all tissue present in the atlas. Data source is from the MethAtlas
#' GitHub repo (\url{https://github.com/nloyfer/meth_atlas}).
#' @usage data(hg19.ref.anc.se)
#' @format An object of class \code{GRanges}.
#'
#' @details All coordinates are in hg19.
#' @references Moss, J., Magenheim, J., Neiman, D. et al. Comprehensive human
#' cell-type methylation atlas reveals origins of circulating cell-free DNA in
#' health and disease. \emph{Nat Commun 9, 5068} (2018).
#' @export
#'
#' @examples
#'   data(hg19.ref.anc.se)
#'   hg19.ref.anc.se
"hg19.ref.anc.se"

#' Hg19 genomic information of cell-type-specific marker CpGs inferred from DNA
#' methylation atlas published by Moss 2018 Nat. Commun.
#' (https://www.nature.com/articles/s41467-018-07466-6).
#'
#'
#' @description Default reference cell type-specific markers used as default in
#' \code{\link{decemedip}}. This dataset represents a GRanges object that contains the collection
#' of Illumina HumanMethylation450K probes that have methylation level less than 0.1 or
#' greater than 0.9 in all tissue present in the atlas. Data source of the methylation atlas is
#' from the MethAtlas GitHub repo (https://github.com/nloyfer/meth_atlas). For details
#' of how the marker CpGs are selected, please refer to the decemedip manuscript.
#' @usage data(hg19.ref.cts.se)
#' @format An object of class \code{GRanges}. rowData(hg19.ref.cts.se) contains
#' information of the selected probes.
#'
#' @details All coordinates are in hg19.
#' @references Moss, J., Magenheim, J., Neiman, D. et al. Comprehensive human
#' cell-type methylation atlas reveals origins of circulating cell-free DNA in
#' health and disease. \emph{Nat Commun 9, 5068} (2018).
#' https://doi.org/10.1038/s41467-018-07466-6
#' @export
#'
#' @examples
#'   data(hg19.ref.cts.se)
#'   hg19.ref.cts.se
"hg19.ref.cts.se"

#' Hg38 genomic information of cell-type-specific marker CpGs inferred from DNA
#' methylation atlas published by Moss 2018 Nat. Commun.
#' (https://www.nature.com/articles/s41467-018-07466-6).
#' @description Same as \code{data(hg19.ref.anc.se)} but lifted over to hg38.
#' @format An object of class \code{GRanges}.
#'
#' @details All coordinates are in hg38.
#' @usage data(hg38.ref.anc.se)
#' @format An object of class \code{GRanges}.
#' @export
#'
#' @examples
#'   data(hg38.ref.anc.se)
#'   hg38.ref.anc.se
"hg38.ref.anc.se"


#' Hg38 genomic information of cell-type-specific marker CpGs inferred from DNA
#' methylation atlas published by Moss 2018 Nat. Commun.
#' (https://www.nature.com/articles/s41467-018-07466-6).
#' @description Same as \code{data(hg19.ref.cts.se)} but lifted over to hg38.
#' @format An object of class \code{GRanges}.
#'
#' @details All coordinates are in hg38.
#' @usage data(hg38.ref.cts.se)
#' @format An object of class \code{GRanges}.
#' @export
#'
#' @examples
#'   data(hg38.ref.cts.se)
#'   hg38.ref.cts.se
"hg38.ref.cts.se"

#' Genomic information of CpG positions in hg19
#'
#' @description This dataset represents a GRanges object that contains the collection
#' of CpGs in chr1-chr22, chrX and chrY. Each row is a CpG. Information was extracted
#' from \code{BSgenome.Hsapiens.UCSC.hg19} of package \code{BSgenome}.
#' @usage data(hg19.cpg.coords)
#' @format An object of class \code{GRanges}.
#' @export
#'
#' @examples
#'   data(hg19.cpg.coords)
#'   hg19.cpg.coords
"hg19.cpg.coords"


#' MeDIP-seq read counts on reference cell type-specific CpGs of 3 PDX samples from Berchuck et al. 2022
#'
#' @description This dataset represents a \code{SummarizedExperiment} object that contains
#' MeDIP-seq read counts on reference cell type-specific CpGs of 3 PDX samples from Berchuck
#' et al. 2022. Each row is a CpG.
#' @usage data(pdx.counts.cts.se)
#' @format An object of class \code{SummarizedExperiment}.
#'
#' @details All coordinates are in hg19.
#' @export
#'
#' @examples
#'   data(pdx.counts.cts.se)
#'   pdx.counts.cts.se
"pdx.counts.cts.se"

#' MeDIP-seq read counts on reference anchor CpGs of 3 PDX samples from Berchuck et al. 2022
#'
#' @description This dataset represents a \code{SummarizedExperiment} object that contains
#' MeDIP-seq read counts on reference anchor CpGs of 3 PDX samples from Berchuck
#' et al. 2022. Each row is a CpG.
#' @usage data(pdx.counts.anc.se)
#' @format An object of class \code{SummarizedExperiment}.
#' @export
#'
#' @details All coordinates are in hg19.
#'
#' @examples
#'   data(pdx.counts.anc.se)
#'   pdx.counts.anc.se
"pdx.counts.anc.se"

