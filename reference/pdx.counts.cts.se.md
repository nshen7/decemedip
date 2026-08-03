# MeDIP-seq read counts on reference cell type-specific CpGs of 3 PDX samples from Berchuck et al. 2022

This dataset represents a `SummarizedExperiment` object that contains
MeDIP-seq read counts on reference cell type-specific CpGs of 3 PDX
samples from Berchuck et al. 2022. Each row is a CpG.

## Usage

``` r
data(pdx.counts.cts.se)
```

## Format

An object of class `SummarizedExperiment`.

## Details

All coordinates are in hg19.

## Examples

``` r
  data(pdx.counts.cts.se)
  pdx.counts.cts.se
#> class: RangedSummarizedExperiment 
#> dim: 2500 3 
#> metadata(0):
#> assays(1): counts
#> rownames(2500): cg18856478 cg20820767 ... cg01071459 cg20726993
#> rowData names(4): probe label pos n_cpgs_100bp
#> colnames(3): LuCaP_147CR LuCap_23.1CR LuCaP_70CR
#> colData names(0):
```
