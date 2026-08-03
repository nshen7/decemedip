# Obtain read counts of regions of interest from MeDIP-seq bam files

Obtain read counts of regions of interest from MeDIP-seq bam files

## Usage

``` r
getRoiReadCount(
  sample_bam_files,
  sample_names,
  sample_paired,
  roi,
  col_data = NULL,
  row_data = NULL,
  bs_genome = "BSgenome.Hsapiens.UCSC.hg19",
  ...
)
```

## Arguments

- sample_bam_files:

  An atomic vector that contains directory paths of the MeDIP-seq sorted
  bam files.

- sample_names:

  An atomic vector of strings that indicates the sample names. Default
  is NULL. If not NULL, please make sure that `sample_names` corresponds
  to elements of `sample_bam_files`.

- sample_paired:

  A logic value that indicates sample_paired-end reads (if TRUE) or
  single-end reads (if FALSE).

- roi:

  A `GRanges` object that contains the genomic coordinates of the region
  of interest (ROI).

- col_data:

  A `DataFrame` object that contains metadata for columns (i.e.,
  samples) if specified. Default is NULL. If not NULL, please make sure
  that rows of `col_data` corresponds to elements of sample_bam_files.
  If input is a non-`DataFrame` object, it will be converted to a
  `DataFrame`.

- row_data:

  A `DataFrame` object that contains metadata for rows (i.e., genomic
  regions) if specified. Default is NULL. If not NULL, please make sure
  that rows of `row_data` corresponds to elements of sample_bam_files.
  If input is a non-`DataFrame` object, it will be converted to a
  `DataFrame`.

- bs_genome:

  A character value that indicates the reference genome name as defined
  by `BSgenome` package. Default is 'BSgenome.Hsapiens.UCSC.hg19'.

- ...:

  Additional arguments passed into
  [`MEDIPS::MEDIPS.createROIset`](https://rdrr.io/pkg/MEDIPS/man/MEDIPS.createROIset.html)

## Value

An object of class `SummarizedExperiment` with read count matrix stored
as an assay named 'counts' (can be extracted using
[`SummarizedExperiment::assays`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html))

## Examples

``` r
# \donttest{
se <- getRoiReadCount(
 sample_bam_files = c('dir/to/bam1','dir/to/bam2'),
 sample_names = c('sample1', 'sample2'),
 sample_paired = TRUE
)
#> Error in getRoiReadCount(sample_bam_files = c("dir/to/bam1", "dir/to/bam2"),     sample_names = c("sample1", "sample2"), sample_paired = TRUE): argument "roi" is missing, with no default
# }
```
