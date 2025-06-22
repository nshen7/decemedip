library(testthat)
library(SummarizedExperiment)

test_that("sample_bam_file and paired_end input are correctly checked", {
  expect_error(decemedip(sample_bam_file = NULL, paired_end = TRUE),
               "Please input 'sample_bam_file'.")

  expect_error(decemedip(sample_bam_file = "sample.bam", paired_end = NULL),
               "Please input 'paired_end'.")
})

test_that("Mutual exclusivity of sample_bam_file and counts inputs", {
  counts_cts <- c(1, 2, 3)
  counts_anc <- c(1, 2, 3)
  expect_error(decemedip(sample_bam_file = "sample.bam", counts_cts = counts_cts,
                         paired_end = TRUE),
               "Invalid: Both counts and bam files are received as input.")
})

test_that("Correctness of counts_cts and counts_anc dimensions", {
  ref_cts <- SummarizedExperiment::SummarizedExperiment(matrix(1:6, ncol=3))
  ref_anc <- SummarizedExperiment::SummarizedExperiment(matrix(1:6, ncol=3))

  expect_error(decemedip(counts_cts = c(1, 2), counts_anc = c(1, 2), ref_cts = ref_cts, ref_anc = ref_anc),
               "Need column 'n_cpgs_100bp' in 'ref_cts' or 'ref_anc', please use 'decemedip::makeReferencePanel' to generate reference panels.")
})

test_that("Integer and numeric checks for counts_cts and counts_anc", {
  ref_cts <- SummarizedExperiment::SummarizedExperiment(matrix(1:6, ncol=3))
  rowData(ref_cts)$n_cpgs_100bp <- 1:2
  ref_anc <- SummarizedExperiment::SummarizedExperiment(matrix(1:6, ncol=3))
  rowData(ref_anc)$n_cpgs_100bp <- 1:2

  expect_error(decemedip(counts_cts = c(1.5, 2.1), counts_anc = c(1, 2), ref_cts = ref_cts, ref_anc = ref_anc),
               "all\\(counts_cts == floor\\(counts_cts\\)\\)")
  expect_error(decemedip(counts_cts = c(1, 2), counts_anc = c(1.5, 2), ref_cts = ref_cts, ref_anc = ref_anc),
               "all\\(counts_anc == floor\\(counts_anc\\)\\)")
})

test_that("ref_assembly argument check", {
  expect_error(decemedip(ref_assembly = 'invalid'),
               "'ref_assembly' can only be hg19 or hg38.")
})

test_that("Reference panel type checks", {
  invalid_obj <- matrix(1:4, ncol = 2)
  expect_error(decemedip(ref_cts = invalid_obj),
               "'ref_cts' and 'ref_anc' have to be RangedSummarizedExperiment or SummarizedExperiment objects.")
})

