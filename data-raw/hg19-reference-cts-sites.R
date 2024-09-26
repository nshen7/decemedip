library(here)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## 450K array probe annotation
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## CpG coordinates in hg19
load('data/hg19.cpg.coords.rda')

## Load full atlas from MethAtlas
### NEED INTERNET
ref_atlas.df <- fread('https://github.com/nloyfer/meth_atlas/raw/master/reference_atlas.csv') |> ## only to obtain colnames
  rename('CpGs' = 'probe')
full_atlas.df <- fread('https://github.com/nloyfer/meth_atlas/raw/master/full_atlas.csv.gz',
                       header = F, col.names = colnames(ref_atlas.df)) |>
  na.omit() |>
  as.data.frame()
full_beta.df <- full_atlas.df[, -1]
rownames(full_beta.df) <- full_atlas.df$probe

## Obtain diff of target and MAX background beta for each cell type
full_diff_max.df <- full_atlas.df |> select('probe')
for (i in seq_len(ncol(full_beta.df)))  {
  label <- colnames(full_beta.df)[i]
  beta_tg <- full_beta.df[[ i]]
  beta_bg <- full_beta.df[, -i] |> as.matrix() |> rowMaxs()
  full_diff_max.df[[label]] <- beta_tg - beta_bg
  print(i)
}

## Obtain diff of target and MEAN background beta for each cell type
full_diff_mean.df <- full_atlas.df |> select('probe')
for (i in seq_len(ncol(full_beta.df)))  {
  label <- colnames(full_beta.df)[i]
  beta_tg <- full_beta.df[[ i]]
  beta_bg <- full_beta.df[, -i] |> as.matrix() |> rowMeans()
  full_diff_mean.df[[label]] <- beta_tg - beta_bg
  print(i)
}

# ---- assemble ref atlas ----
K_PER_CT <- 100
K_ADD <- 500; N_PER_SET <- 10

getCTSHyperProbes <- function(label, top_k) {
  idx <- which(full_diff_max.df[[label]] > 0)
  df <- data.frame(probe = full_atlas.df[['probe']][idx],
                   label = paste0(label, ' hypermethylated'),
                   diff_mean = full_diff_mean.df[[label]][idx],
                   diff_max  = full_diff_max.df[[label]][idx]) |>
    group_by(label) |>
    slice_max(order_by = diff_mean, n = top_k)
  return(df)
}

## Select cell type specific marker probes that are top-ranked in diff between
## target cell type and average background
cts_probes.temp <- map_dfr(colnames(full_diff_max.df)[-1],
                           getCTSHyperProbes, top_k = K_PER_CT) |>
  select(probe, label)
ref_cts_beta.temp <- full_beta.df[match(cts_probes.temp$probe, full_atlas.df$probe), ]

## Iteratively identifying the two most similar cell types in the atlas,
## and adding the CpG site upon which these two cell types differ the most
left_atlas.temp <- full_atlas.df |> filter(!probe %in% cts_probes.temp$probe)

for (i in 1:(K_ADD/N_PER_SET)) {
  dist.df <- dist(t(ref_cts_beta.temp), method = 'manhattan') |>
    as.matrix() |>
    as.data.frame() |>
    mutate(cell_type_1 = colnames(ref_cts_beta.temp)) |>
    pivot_longer(cols = -c('cell_type_1'),
                 values_to = 'distance',
                 names_to = 'cell_type_2') |>
    filter(distance > 0)

  # Find the two cell types with lowest distance
  idx <- which.min(dist.df$distance)
  c1 <- dist.df$cell_type_1[idx]
  c2 <- dist.df$cell_type_2[idx]

  # Find 10 probes that distinguish them, 5 from each direction
  probes_to_add_1 <- left_atlas.temp |>
    mutate(diff = left_atlas.temp[c1] - left_atlas.temp[c2]) |>
    slice_max(order_by = diff, n = N_PER_SET/2, with_ties = FALSE) |>
    select(probe) |>
    mutate(label = 'Additional set')
  probes_to_add_2 <- left_atlas.temp |>
    mutate(diff = left_atlas.temp[c2] - left_atlas.temp[c1]) |>
    slice_max(order_by = diff, n = N_PER_SET/2, with_ties = FALSE) |>
    select(probe) |>
    mutate(label = 'Additional set')
  probes_to_add <- rbind(probes_to_add_1, probes_to_add_2)

  # Add the probes to CTS reference and remove them from left_atlas.temp
  cts_probes.temp <- rbind(cts_probes.temp, probes_to_add)
  ref_cts_beta.temp <- rbind(ref_cts_beta.temp,
                             full_beta.df[match(probes_to_add$probe, full_atlas.df$probe), ])
  left_atlas.temp <- full_atlas.df |> filter(!probe %in% cts_probes.temp$probe)

}

## The complete reference panel
ref_cts.df <- cbind(
  cts_probes.temp,
  as.data.frame(anno450k[match(cts_probes.temp$probe, rownames(anno450k)), c("chr", "pos")])
) |>
  mutate(start = pos, end = pos) |>
  as.data.frame()
rownames(ref_cts.df) <- ref_cts.df$probe

ref_cts_beta.df <- full_beta.df[match(cts_probes.temp$probe, rownames(full_beta.df)), ]

## Sanity check
sum(duplicated(ref_cts.df)) # = 0; no duplicated probes
all(rownames(ref_cts.df) == ref_cts.df$probe) # = TRUE; all probes correctly mapped
all(rownames(ref_cts_beta.df) == ref_cts.df$probe) # = TRUE; all probes correctly mapped


## Write out
ref_cts.gr <- makeGRangesFromDataFrame(ref_cts.df, keep.extra.columns = T)
ref_cts.gr$n_cpgs_100bp <- countOverlaps(ref_cts.gr |> resize(width = 100, fix = 'center'), hg19.cpg.coords)
ref_cts.se <- SummarizedExperiment(rowData = ref_cts.gr,
                                   assays = list(beta = ref_cts_beta.df))
genome(ref_cts.se) <- "hg19"

hg19.ref.cts.se <- ref_cts.se
usethis::use_data(hg19.ref.cts.se, overwrite = TRUE)
usethis::use_data(hg19.ref.cts.se, overwrite = TRUE, internal = TRUE)

