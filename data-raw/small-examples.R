library(here)

## CpG coordinates and reference sites in hg19
load("data/hg19.ref.cts.se.rda")
load("data/hg19.ref.anc.se.rda")
load("data/pdx.counts.cts.se.rda")
load("data/pdx.counts.anc.se.rda")

cell_types <- c(1:7, 22)

## Extract the CpG signatures for blood cell types and prostate only
idx_cts <- c(1:700, 2201:2300)
example.hg19.ref.cts.se <- hg19.ref.cts.se[idx_cts, cell_types]
example.pdx.counts.cts.se <- pdx.counts.cts.se[idx_cts, ]

## Reduce the set of anchor CpGs
idx_anc <- c(1:150, 501:650)
example.hg19.ref.anc.se <- hg19.ref.anc.se[idx_anc, cell_types]
example.pdx.counts.anc.se <- pdx.counts.anc.se[idx_anc, ]

usethis::use_data(example.hg19.ref.cts.se, overwrite = TRUE)
usethis::use_data(example.pdx.counts.cts.se, overwrite = TRUE)

usethis::use_data(example.hg19.ref.anc.se, overwrite = TRUE)
usethis::use_data(example.pdx.counts.anc.se, overwrite = TRUE)



