#' Main function to perform cell type deconvolution with MeDIP-seq data
#'
#' @param sample_bam_file A string value that indicates the file path to the bam file of a
#' MeDIP-seq sample of interest. If `sample_bam_file` is specified, please do not
#' specify `counts_cts` and `counts_anc` to avoid conflict.
#' @param paired_end A logic value that indicates whether the bam file is from paired-end
#' reads or single-end. Specify `TRUE} for paired-end and \code{FALSE` for single-end.
#' @param counts_cts An atomic vector of integer values that indicates the read counts of
#' a MeDIP-seq sample on reference sites/regions. If `counts_cts` and `counts_anc`
#' is specified, please do not specify `sample_bam_file` and `paired_end` to avoid conflict.
#' @param counts_anc An atomic vector of integer values that indicates the read counts of
#' a MeDIP-seq sample on reference sites/regions. If  `counts_cts` and `counts_anc`
#' is specified, please do not specify `sample_bam_file` and `paired_end` to avoid conflict.
#' @param diagnostics A logic value that indicates whether to include components of the stan model
#' in the output that are necessary for future diagnostics of the model, such as posterior
#' predictive checks. For details, please refer to the function `\link{plotDiagnostics}`.
#' @param ref_assembly A string that represents the genome assembly that should be used for cell type-specific
#' sites in the reference panel ('hg19' or 'hg38'). Default to 'hg19'. The default reference is explained
#' in the manuscript of \pkg{decemedip}. Alternatively, if the user want to provide their own
#' reference panel by using the `ref_cts` and `ref_anc` arguments.
#' @param ref_cts A `SummarizedExperiment` object that contains the genomic coordinates and
#' beta values of the cell type-specific sites/regions from reference cell types.
#' The `\link{makeReferencePanel}` can be used to generate such a panel.
#' @param ref_anc Same as `ref_cts` but for anchor sites.
#' @param weight_cts A numeric value indicating the weights that should be put on cell type-specific
#' sites/regions. Default is 0.5.
#' @param weight_anc A numeric value indicating the weights that should be put on cell type-specific
#' sites/regions. Default is 1.
#' @param seed The seed for random number generation in MCMC sampling.
#' @param chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param iter A positive integer specifying the number of iterations for each chain (including
#' warmup). The default is 2000.
#' @param cores A positive integer specifying the number of cores that can be used for
#' MCMC sampling. The default is 1.
#' @param stan_input_params A named list of parameters that specifies the prior of the
#' decemedip model.
#' @param stan_control A named list of parameters to control the sampler's behavior in
#' Stan. See the details in the documentation for the control argument in `\link[rstan]{stan}`.
#' @param timeout_sec A numerical value indicating the CPU/processor time (in seconds) allowed for
#' the MCMC to run before restarting the chains with a new random seed.
#' @param max_retries An integer value indicating the maximum number of tries with different seed
#' for MCMC if it fails to converge.
#' @param ... Other parameters that can be passed to the `\link[rstan]{sampling}` function.
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment ncol
#' @importFrom SummarizedExperiment nrow
#' @importFrom R.utils withTimeout
#' @return A list of two elements:
#' 1. `data_list`: An organized list of variables used as input to the Stan posterior sampling function.
#' 2. `posterior`: An `stanfit` object produced by Stan representing the fitted posteriors.
#' @export
#'
#' @examples
#' data(pdx.counts.cts.se)
#' data(pdx.counts.anc.se)
#' # read counts of cell type-specific CpGs of the sample 'LuCaP_147CR'
#' counts_cts <- SummarizedExperiment::assays(pdx.counts.cts.se)$counts[, "LuCaP_147CR"]
#' # read counts of anchor CpGs of the sample 'LuCaP_147CR'
#' counts_anc <- SummarizedExperiment::assays(pdx.counts.anc.se)$counts[, "LuCaP_147CR"]
#' \dontrun{
#' # Fit decemedip model (iter=100 for demonstration, by default iter=2000)
#' \donttest{
#' output <- decemedip(counts_cts = counts_cts, counts_anc = counts_anc, iter = 100, cores = 1, chains = 1)
#' }
#'
decemedip <- function(
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
    cores = 1,
    chains = 4,
    iter = 2000,
    stan_input_params = list(
      "s_mu"      = 3,
      "s_sigma"   = 3,
      "n_knot_z"  = 0,
      "degree_z"  = 3,
      "Xi"        = cor(as.matrix(SummarizedExperiment::assay((ref_cts)))),
      "s_theta"   = 3,
      "s_tau"     = 3
    ),
    stan_control = list(adapt_delta = 0.95, max_treedepth = 15),
    timeout_sec = 2 * (diagnostics + 1) * iter * chains,
    max_retries = 3,
    ...) {
  temp_env <- new.env()

  if (is.null(ref_cts) & is.null(ref_anc)) {
    if (ref_assembly == "hg19") {
      data(hg19.ref.cts.se, envir = temp_env)
      ref_cts <- temp_env$hg19.ref.cts.se
      data(hg19.ref.anc.se, envir = temp_env)
      ref_anc <- temp_env$hg19.ref.anc.se
    } else if (ref_assembly == "hg38") {
      data(hg38.ref.cts.se, envir = temp_env)
      ref_cts <- temp_env$hg38.ref.cts.se
      data(hg38.ref.anc.se, envir = temp_env)
      ref_anc <- temp_env$hg38.ref.anc.se
    } else {
      stop("'ref_assembly' can only be hg19 or hg38.")
    }
  }

  ## Checks on reference
  if ((!is(ref_cts, "RangedSummarizedExperiment") & !is(ref_cts, "SummarizedExperiment")) |
    (!is(ref_anc, "RangedSummarizedExperiment") & !is(ref_anc, "SummarizedExperiment"))) {
    stop("'ref_cts' and 'ref_anc' have to be RangedSummarizedExperiment or SummarizedExperiment objects.")
  }

  if ((!"n_cpgs_100bp" %in% colnames(rowData(ref_cts))) | (!"n_cpgs_100bp" %in% colnames(rowData(ref_anc)))) {
    stop("Need column 'n_cpgs_100bp' in 'ref_cts' or 'ref_anc', please use 'decemedip::makeReferencePanel' to generate reference panels.")
  }

  ## Checks on stan parameters
  stopifnot(c("s_mu", "s_sigma", "n_knot_z", "degree_z", "Xi", "s_theta", "s_tau") %in% names(stan_input_params))

  ## Checks on exclusivity of counts and bam files
  if (!is.null(sample_bam_file) | !is.null(paired_end)) {
    if (is.null(sample_bam_file)) stop("Please input 'sample_bam_file'")
    if (is.null(paired_end)) stop("Please input 'paired_end'")
    if (length(counts_cts) != 0) stop("Invalid: Both counts and bam files are received as input.")
    if (length(counts_anc) != 0) stop("Invalid: Both counts and bam files are received as input.")
  }
  if (length(counts_cts) > 0 | length(counts_anc) > 0) {
    # check whether dim aligns with reference
    stopifnot(length(counts_cts) == nrow(ref_cts))
    stopifnot(length(counts_anc) == nrow(ref_anc))
    # check whether values are integers
    stopifnot(all(counts_cts == floor(counts_cts)) && is.numeric(counts_cts))
    stopifnot(all(counts_anc == floor(counts_anc)) && is.numeric(counts_anc))
    # check on exclusivity of counts and bam files
    if (!is.null(sample_bam_file)) stop("Invalid: Both counts and bam files are received as input.")
    if (!is.null(paired_end)) stop("Invalid: Both counts and bam files are received as input.")
  }

  ## Load read counts of reference sites if bam file is provided
  if (!is.null(sample_bam_file) & !is.null(paired_end)) {
    sample_cts.se <- getRoiReadCount(sample_bam_file, sample_names = NULL, sample_paired = paired_end, roi = granges(ref_cts))
    sample_anc.se <- getRoiReadCount(sample_bam_file, sample_names = NULL, sample_paired = paired_end, roi = granges(ref_anc))
    counts_cts <- SummarizedExperiment::assay((sample_cts.se))[, 1] |> unlist()
    counts_anc <- SummarizedExperiment::assay((sample_anc.se))[, 1] |> unlist()
  }

  ## Check for low coverage or unusual high counts
  if (sum(counts_cts) < 0.25 * length(counts_cts)) {
    warning("Unusual low coverage for cell type-specific reference sites.", immediate. = TRUE)
  }

  ## Prepare model variables
  y <- c(counts_cts, counts_anc)
  X <- rbind(SummarizedExperiment::assay(ref_cts), SummarizedExperiment::assay((ref_anc)))
  z <- c(SummarizedExperiment::rowData(ref_cts)$n_cpgs_100bp, SummarizedExperiment::rowData(ref_anc)$n_cpgs_100bp) |> log1p()
  ## Assign weights to the CTS sites and anchor sites
  weights_raw <- c(rep(weight_cts, SummarizedExperiment::nrow(ref_cts)), rep(weight_anc, SummarizedExperiment::nrow(ref_anc)))

  # add head and tail into knots for computation in the stan script
  if (stan_input_params$n_knot_z == 0) probs <- NULL else probs <- seq(0, 1, length.out = stan_input_params$n_knot_z)
  stan_input_params$knots_z <- quantile(z, probs = c(0, probs, 1))
  stan_input_params$n_knot_z <- length(stan_input_params$knots_z)

  ## Prepare input for stan model
  data_list <- c(
    list("N" = nrow(X), "K" = ncol(X), "y" = y, "X" = X, "z" = z, "weights_raw" = weights_raw),
    stan_input_params
  )


  if (diagnostics) {
    model <- stanmodels$decemedip1
  } else {
    ## Run the stan model that does not generate y_sim
    model <- stanmodels$decemedip0
  }

  # if (is.null(stan_control)) stan_control <- list(adapt_delta = 0.95, max_treedepth = 15)

  ## Run MCMC (automatically rerun with a new seed if reaches timeout)
  success <- FALSE # Flag to track success

  for (i in seq_len(max_retries)) {
    this_seed <- seed + (i - 1) * exp(10)
    posterior <- tryCatch(
      {
        R.utils::withTimeout(
          {
            rstan::sampling(model,
              data = data_list,
              iter = iter, chains = chains, cores = cores, seed = this_seed,
              control = stan_control, ...
            )
          },
          timeout = timeout_sec
        ) # Timeout in seconds
      },
      TimeoutException = function(e) {
        message(paste("MCMC with seed", this_seed, "took too long, trying a new one..."))
        NULL # Return NULL if timeout occurs
      }
    )

    if (!is.null(posterior)) {
      success <- TRUE
      message(paste("MCMC converged with seed", this_seed))
      break # Exit the loop if the algorithm succeeds
    }
  }

  if (!success) {
    message(paste("MCMC failed to converge after", max_retries, "tries."))
  }

  ## Return model
  return(list("data_list" = data_list, "posterior" = posterior))
}
