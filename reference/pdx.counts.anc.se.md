# MeDIP-seq read counts on reference anchor CpGs of 3 PDX samples from Berchuck et al. 2022

This dataset represents a `SummarizedExperiment` object that contains
MeDIP-seq read counts on reference anchor CpGs of 3 PDX samples from
Berchuck et al. 2022. Each row is a CpG.

## Usage

``` r
data(pdx.counts.anc.se)
```

## Format

An object of class `SummarizedExperiment`.

## Details

All coordinates are in hg19.

## Examples

``` r
  data(pdx.counts.anc.se)
  pdx.counts.anc.se
#> class: RangedSummarizedExperiment 
#> dim: 1000 3 
#> metadata(0):
#> assays(1): counts
#> rownames(1000): 353534 76294 ... 73948 87963
#> rowData names(7): probe pos ... avg_beta_rank n_cpgs_100bp
#> colnames(3): LuCaP_147CR LuCap_23.1CR LuCaP_70CR
#> colData names(0):
```
