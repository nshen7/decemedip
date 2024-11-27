library(here)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## 450K array probe annotation
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## CpG coordinates and CTS reference sites in hg19
load('data/hg19.cpg.coords.rda')
load('data/hg19.ref.cts.se.rda')

ref_cts.gr <- granges(hg19.ref.cts.se)

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

# ---- Find anchor sites (i.e., all-tissue M/U sites) from MethAtlas ----

K_PER_DIRECTION <- 500

# Find blocks where WGBS samples are all methylated (beta > 1-margin) or all unmethylated (beta < margin)
margin <- 0.1
idx_all_m <- which(rowMeans(full_beta.df > 1-margin) >= 1)
idx_all_u <- which(rowMeans(full_beta.df < margin) >= 1)

buildSE <- function(idx, label) {

  subset.df <- full_atlas.df[idx, ]
  beta.df   <- subset.df[, -1]

  ## Compute average beta values per marker
  avg_beta  <- rowMeans(beta.df)
  if (label == 'All-tissue U') avg_beta_rank <- rank(avg_beta) else avg_beta_rank <- rank(-avg_beta)

  probe_list <- unlist(subset.df[, 1])
  coords     <- anno450k[match(probe_list, rownames(anno450k)), c("chr", "pos")]
  md.df      <- data.frame(probe = probe_list, coords) %>%
    mutate(label = label, margin = margin,
           start = pos, end = pos,
           avg_beta = avg_beta,
           avg_beta_rank = avg_beta_rank)
  rownames(md.df) <- NULL

  se <- SummarizedExperiment(rowData = makeGRangesFromDataFrame(md.df, keep.extra.columns = TRUE),
                             assays = list(beta = beta.df))
  return(se)
}

u_atlas.se <- buildSE(idx = idx_all_u, label = 'All-tissue U')
rowData(u_atlas.se)$n_cpgs_100bp <- countOverlaps(granges(u_atlas.se) |> resize(width = 100, fix = 'center'), hg19.cpg.coords)
m_atlas.se <- buildSE(idx = idx_all_m, label = 'All-tissue M')
rowData(m_atlas.se)$n_cpgs_100bp <- countOverlaps(granges(m_atlas.se) |> resize(width = 100, fix = 'center'), hg19.cpg.coords)

table(rowData(u_atlas.se)$n_cpgs_100bp)
table(rowData(m_atlas.se)$n_cpgs_100bp)

## Sample the anchor probes to follow similar distribution of CpG density as the CTS probes
target_freqs_0 <- data.frame(table(ref_cts.gr$n_cpgs_100bp) / length(ref_cts.gr)) |>
  rename('Var1' = 'n_cpgs_100bp', 'Freq' = 'prob') |>
  mutate(n_cpgs_100bp = as.integer(n_cpgs_100bp))

sampleRowsAccordingToCGDensity <- function(se, k, seed = 2022) {

  set.seed(seed)

  # Freqs in u_atlas.se/m_atlas.se region set
  original_freqs <- data.frame(table(rowData(se)$n_cpgs_100bp)) |>
    rename('Var1' = 'n_cpgs_100bp', 'Freq' = 'count') |>
    mutate(n_cpgs_100bp = as.integer(n_cpgs_100bp))

  target_freqs <- target_freqs_0 |>
    left_join(original_freqs, by = 'n_cpgs_100bp') |>
    mutate(adjusted_prob = prob / count)
  sampling_prob <- rowData(se) |> as.data.frame() |> left_join(target_freqs, by = 'n_cpgs_100bp') |> pull(adjusted_prob)
  sampling_prob[is.na(sampling_prob)] <- 0
  sampled.se <- se[sample(x = 1:nrow(se), size = k, prob = sampling_prob, replace = FALSE), ]

  return(sampled.se)
}

u_atlas_sampled.se <- sampleRowsAccordingToCGDensity(se = u_atlas.se, k = K_PER_DIRECTION)
m_atlas_sampled.se <- sampleRowsAccordingToCGDensity(se = m_atlas.se, k = K_PER_DIRECTION)

table(rowData(u_atlas_sampled.se)$n_cpgs_100bp)
table(rowData(m_atlas_sampled.se)$n_cpgs_100bp)


ref_anc.se <- rbind(u_atlas_sampled.se, m_atlas_sampled.se)

genome(ref_anc.se) <- "hg19"

hg19.ref.anc.se <- ref_anc.se
usethis::use_data(hg19.ref.anc.se, overwrite = TRUE)
usethis::use_data(hg19.ref.anc.se, overwrite = TRUE, internal = TRUE)
