library(here)
library(tidyverse)
library(data.table)
library(decemedip)
devtools::load_all() # to load new changes in functions

md_samples <- fread(here('..', 'cfMeDIP-deconv-experiments', 'data', 'metadata',
                         'Berchuck2022_LuCap_PDX_MeDIP',
                         'sample_metadata_processed_Berchuck2022_LuCap_PDX_MeDIP.csv'))
idx <- 1:3
pdx.counts.cts.se <- getRoiReadCount(sample_bam_files = md_samples$bam_dir[idx],
                                 sample_names = md_samples$Sample_Name[idx],
                                 sample_paired = TRUE,
                                 roi = granges(hg19.ref.cts.se))
pdx.counts.anc.se <- getRoiReadCount(sample_bam_files = md_samples$bam_dir[idx],
                                 sample_names = md_samples$Sample_Name[idx],
                                 sample_paired = TRUE,
                                 roi = granges(hg19.ref.anc.se))
usethis::use_data(pdx.counts.cts.se, pdx.counts.anc.se, overwrite = TRUE)
