# Main function to perform cell type deconvolution with MeDIP-seq data

Main function to perform cell type deconvolution with MeDIP-seq data

## Usage

``` r
decemedip(
  sample_bam_file = NULL,
  paired_end = NULL,
  counts_cts = c(),
  counts_anc = c(),
  ref_assembly = "hg19",
  ref_cts = NULL,
  ref_anc = NULL,
  weight_cts = 1,
  weight_anc = 0.5,
  diagnostics = TRUE,
  seed = 2024,
  cores = 4,
  chains = 4,
  iter = 2000,
  stan_input_params = list(s_mu = 3, s_sigma = 3, n_knot_z = 0, degree_z = 3, Xi =
    cor(as.matrix(SummarizedExperiment::assays(ref_cts)[[1]])), s_theta = 3, s_tau = 3),
  stan_control = list(adapt_delta = 0.95, max_treedepth = 15),
  timeout_sec = 2 * (diagnostics + 1) * iter * chains,
  max_retries = 3,
  ...
)
```

## Arguments

- sample_bam_file:

  A string value that indicates the file path to the bam file of a
  MeDIP-seq sample of interest. If `sample_bam_file` is specified,
  please do not specify `counts_cts` and `counts_anc` to avoid conflict.

- paired_end:

  A logic value that indicates whether the bam file is from paired-end
  reads or single-end. Specify `TRUE` for paired-end and `FALSE` for
  single-end.

- counts_cts:

  An atomic vector of integer values that indicates the read counts of a
  MeDIP-seq sample on reference sites/regions. If `counts_cts` and
  `counts_anc` is specified, please do not specify `sample_bam_file` and
  `paired_end` to avoid conflict.

- counts_anc:

  An atomic vector of integer values that indicates the read counts of a
  MeDIP-seq sample on reference sites/regions. If `counts_cts` and
  `counts_anc` is specified, please do not specify `sample_bam_file` and
  `paired_end` to avoid conflict.

- ref_assembly:

  A string that represents the genome assembly that should be used for
  cell type-specific sites in the reference panel ('hg19' or 'hg38').
  Default to 'hg19'. The default reference is explained in the
  manuscript of decemedip. Alternatively, if the user want to provide
  their own reference panel by using the `ref_cts` and `ref_anc`
  arguments.

- ref_cts:

  A `SummarizedExperiment` object that contains the genomic coordinates
  and beta values of the cell type-specific sites/regions from reference
  cell types. The
  [`makeReferencePanel`](https://nshen7.github.io/decemedip/reference/makeReferencePanel.md)
  can be used to generate such a panel.

- ref_anc:

  Same as `ref_cts` but for anchor sites.

- weight_cts:

  A numeric value indicating the weights that should be put on cell
  type-specific sites/regions. Default is 0.5.

- weight_anc:

  A numeric value indicating the weights that should be put on cell
  type-specific sites/regions. Default is 1.

- diagnostics:

  A logic value that indicates whether to include components of the stan
  model in the output that are necessary for future diagnostics of the
  model, such as posterior predictive checks. For details, please refer
  to the function
  [`plotDiagnostics`](https://nshen7.github.io/decemedip/reference/plotDiagnostics.md).

- seed:

  The seed for random number generation in MCMC sampling.

- cores:

  A positive integer specifying the number of cores that can be used for
  MCMC sampling. The default is 4

- chains:

  A positive integer specifying the number of Markov chains. The default
  is 4.

- iter:

  A positive integer specifying the number of iterations for each chain
  (including warmup). The default is 2000.

- stan_input_params:

  A named list of parameters that specifies the prior of the decemedip
  model.

- stan_control:

  A named list of parameters to control the sampler's behavior in Stan.
  See the details in the documentation for the control argument in
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- timeout_sec:

  A numerical value indicating the CPU/processor time (in seconds)
  allowed for the MCMC to run before restarting the chains with a new
  random seed.

- max_retries:

  An integer value indicating the maximum number of tries with different
  seed for MCMC if it fails to converge.

- ...:

  Other parameters that can be passed to the
  [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  function.

## Value

A list of two elements:

1.  `data_list`: An organized list of variables used as input to the
    Stan posterior sampling function.

2.  `posterior`: An `stanfit` object produced by Stan representing the
    fitted posteriors.

## Examples

``` r
data(pdx.counts.cts.se)
data(pdx.counts.anc.se)
# read counts of cell type-specific CpGs of the sample 'LuCaP_147CR'
counts_cts <- SummarizedExperiment::assays(pdx.counts.cts.se)$counts[,'LuCaP_147CR']
# read counts of anchor CpGs of the sample 'LuCaP_147CR'
counts_anc <- SummarizedExperiment::assays(pdx.counts.anc.se)$counts[,'LuCaP_147CR']
# Fit decemedip model (iter=100 for demonstration, by default iter=2000)
output <- decemedip(counts_cts = counts_cts, counts_anc = counts_anc, iter = 100, cores = 1, chains = 1)
#> 
#> SAMPLING FOR MODEL 'decemedip1' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.002351 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 23.51 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 7
#> Chain 1:            adapt_window = 38
#> Chain 1:            term_buffer = 5
#> Chain 1: 
#> Chain 1: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 1: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 1: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 1: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 1: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 1: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 1: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 1: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 1: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 1: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 1: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 1: Iteration: 100 / 100 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 33.034 seconds (Warm-up)
#> Chain 1:                21.824 seconds (Sampling)
#> Chain 1:                54.858 seconds (Total)
#> Chain 1: 
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> MCMC converged with seed 2024
```
