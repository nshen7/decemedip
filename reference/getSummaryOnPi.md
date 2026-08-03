# Extract summary statistics and diagnostics on fitted cell type proportions

Extract summary statistics and diagnostics on fitted cell type
proportions

## Usage

``` r
getSummaryOnPi(
  posterior,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  digits_summary = 5,
  cell_type_names = NULL,
  ...
)
```

## Arguments

- posterior:

  The fitted posterior object from decemedip model.

- probs:

  A numeric vector specifying quantiles of interest. The defaults is
  c(0.025,0.25,0.5,0.75,0.975).

- digits_summary:

  The number of significant digits to use in the summary, defaulting to
  5.

- cell_type_names:

  Name of the cell types in reference panel. The order should align with
  order in the reference panel.

- ...:

  Additional arguments that get passed into
  [`rstan::monitor`](https://mc-stan.org/rstan/reference/monitor.html)
  function.

## Value

A data.frame object containg summary statistics and diagnostic
statistics of the fitted cell type proportions.

## Examples

``` r

data(pdx.counts.cts.se)
data(pdx.counts.anc.se)
# read counts of cell type-specific CpGs of the sample 'LuCaP_147CR'
counts_cts <- SummarizedExperiment::assays(pdx.counts.cts.se)$counts[,'LuCaP_147CR']
# read counts of anchor CpGs of the sample 'LuCaP_147CR'
counts_anc <- SummarizedExperiment::assays(pdx.counts.anc.se)$counts[,'LuCaP_147CR']
# Fit decemedip model (iter=100 for demonstration, by default iter=2000)
output <- decemedip(counts_cts = counts_cts, counts_anc = counts_anc, iter = 100)
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

smr_pi.df <- getSummaryOnPi(output$posterior)
```
