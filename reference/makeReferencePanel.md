# Function for assembling a SummarizedExperiment object of reference panel in **hg19** for cell type deconvolution (which is used in `decemedip` function)

Function for assembling a SummarizedExperiment object of reference panel
in **hg19** for cell type deconvolution (which is used in `decemedip`
function)

## Usage

``` r
makeReferencePanel(
  row_ranges,
  X,
  cpg_coords,
  col_names = NULL,
  row_names = NULL,
  col_data = NULL,
  row_data = NULL
)
```

## Arguments

- row_ranges:

  A `GRanges` object that contains the genomic coordinates of *reference
  regions/sites*.

- X:

  A matrix that contains the beta values of reference regions. **Each
  row is a region and each column is a cell type.** If `X` has row names
  or column names, the output `SummarizedExperiment` object will inherit
  them.

- cpg_coords:

  A `GRanges` object that contains genomic coordinates of all CpGs in
  the genome. A ready-to-use CpG list for hg19 is available to download
  at
  <https://github.com/nshen7/decemedip-experiments/blob/main/hg19.cpg.coords.rda>.
  It is used for generating coloum `n_cpg_100bp` in the reference panel,
  which represents CpG density around the reference CpG.

- col_names:

  An atomic vector of strings that indicates the column names, i.e.,
  names of the cell types. Default is NULL. If not NULL, the column
  names of `X` will be overwritten by this argument.

- row_names:

  An atomic vector of strings that indicates the row names, i.e., names
  of the reference regions. Default is NULL. If not NULL, the row names
  of `X` will be overwritten by this argument.

- col_data:

  A `DataFrame` object that contains metadata for columns (i.e., cell
  types) if specified. Each row in `col_data` should contain info of a
  cell type in `X`. If input is a non-`DataFrame` object, it will be
  converted to a `DataFrame`. Default is NULL.

- row_data:

  A `DataFrame` object that contains metadata for row (i.e., reference
  regions) if specified. Each row in `row_data` should contain info of a
  reference region in `X`. If input is a non-`DataFrame` object, it will
  be converted to a `DataFrame`. Default is NULL.

## Value

An `SummarizedExperiment` object with each row represents a reference
region and an assay named 'X' that stores the beta values of reference
regions.

## Examples

``` r

row_ranges <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr3")),
  ranges = IRanges::IRanges(start = c(100, 200, 300),
                            end = c(100, 200, 300)),
  cpg_id = c("cpg_1", "cpg_2", "cpg_3")      # CpG site IDs
)

cpg_coords = GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3")),
  ranges = IRanges::IRanges(start = c(100, 101, 200, 201, 300, 301),
                            end = c(100, 101, 200, 201, 300, 301))
)

X = matrix(runif(6), nrow = 3)

makeReferencePanel(
  row_ranges = row_ranges,
  X = X,
  cpg_coords = cpg_coords
)
#> class: RangedSummarizedExperiment 
#> dim: 3 2 
#> metadata(0):
#> assays(1): beta
#> rownames: NULL
#> rowData names(2): cpg_id n_cpgs_100bp
#> colnames: NULL
#> colData names(0):
```
