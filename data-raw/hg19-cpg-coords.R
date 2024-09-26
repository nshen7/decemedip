library(BSgenome)
library(GenomicRanges)
library(tidyverse)
library(data.table)

CHR_NUMS <- as.character(c(1:22, 'X', 'Y'))
CHR_NAMES <- paste0('chr', CHR_NUMS)

# CpG coordinates in reference genome
genome <- BSgenome.Hsapiens.UCSC.hg19
hg19.cpg.coords <- map_dfr(
  CHR_NAMES,
  ~ matchPattern('CG', genome[[.x]]) %>% as.data.frame %>% mutate(chr = .x) %>% select(chr, start, end)
)
hg19.cpg.coords <- GenomicRanges::makeGRangesFromDataFrame(hg19.cpg.coords)

usethis::use_data(hg19.cpg.coords, overwrite = TRUE)
tools::resaveRdaFiles('data/hg19.cpg.coords.rda')
