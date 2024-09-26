library(here)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(rtracklayer)

## 450K array probe annotation
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## CpG coordinates and reference sites in hg19
load('data/hg19.cpg.coords.rda')
load('data/hg19.ref.cts.se.rda')
load('data/hg19.ref.anc.se.rda')

# # Download (run once) and import the chain file
# URL <- 'https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz'
# download.file(URL, destfile = 'data-raw/hg19ToHg38.over.chain.gz', method="curl")
# R.utils::gunzip('data-raw/hg19ToHg38.over.chain.gz'), remove = TRUE, overwrite = TRUE)
chain_file <- 'data-raw/hg19ToHg38.over.chain'
chain <- import.chain(chain_file)

# Lift over the reference probes
ref_cts_lifted.list <- liftOver(granges(hg19.ref.cts.se), chain)
idx_unmapped_cts <- which(elementNROWS(ref_cts_lifted.list) == 0)

ref_anc_lifted.list <- liftOver(granges(hg19.ref.anc.se), chain)
idx_unmapped_anc <- which(elementNROWS(ref_anc_lifted.list) == 0)

# Remove unmapped probes
ref_cts_2.se <- hg19.ref.cts.se; ref_anc_2.se <- hg19.ref.anc.se
if (length(idx_unmapped_cts) > 0) ref_cts_2.se <- ref_cts_2.se[-idx_unmapped_cts, ]
if (length(idx_unmapped_anc) > 0) ref_anc_2.se <- ref_anc_2.se[-idx_unmapped_ans, ]

# Output lifted reference panel
ref_cts_lifted.se <- SummarizedExperiment(
  rowData = unlist(ref_cts_lifted.list),
  assays = list(beta = assays(ref_cts_2.se)$beta)
)
genome(ref_cts_lifted.se) <- "hg38"

hg38.ref.cts.se <- ref_cts_lifted.se
usethis::use_data(hg38.ref.cts.se, overwrite = TRUE)
usethis::use_data(hg38.ref.cts.se, overwrite = TRUE, internal = TRUE)

ref_anc_lifted.se <- SummarizedExperiment(
  rowData = unlist(ref_anc_lifted.list),
  assays = list(beta = assays(ref_anc_2.se)$beta)
)
genome(ref_anc_lifted.se) <- "hg38"

hg38.ref.anc.se <- ref_anc_lifted.se
usethis::use_data(hg38.ref.anc.se, overwrite = TRUE)
usethis::use_data(hg38.ref.anc.se, overwrite = TRUE, internal = TRUE)

