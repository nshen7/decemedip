#' Perform cell type deconvolution with MeDIP-seq data
#'
#' @param sample_bam_file A string value that indicates the file path to the bam file of a
#' MeDIP-seq sample of interest. If \code{sample_bam_file} is specified, please do not
#' specify \code{counts_cts} and \code{counts_anc} to avoid conflict.
#' @param paired_end A logic value that indicates whether the bam file is from paired-end
#' reads or single-end. Specify \code{TRUE} for paired-end and \code{FALSE} for single-end.
#' @param counts_cts An atomic vector of integer values that indicates the read counts of
#' a MeDIP-seq sample on reference sites/regions. If \code{counts_cts} and \code{counts_anc}
#' is specified, please do not specify \code{sample_bam_file} and \code{paired_end} to avoid conflict.
#' @param counts_anc An atomic vector of integer values that indicates the read counts of
#' a MeDIP-seq sample on reference sites/regions. If  \code{counts_cts} and \code{counts_anc}
#' is specified, please do not specify \code{sample_bam_file} and \code{paired_end} to avoid conflict.
#' @param diagnostics A logic value that indicates whether to include components of the stan model
#' in the output that are necessary for future diagnostics of the model, such as posterior
#' predictive checks. For details, please refer to the function \code{\link{plotDiagnostics}}.
#' @param ref_cts A \code{SummarizedExperiment} object that contains the genomic coordinates
#' and beta values of the cell type-specific sites/regions from reference cell types. The default
#' reference is explained in the manuscript of \pkg{decemedip}.
#' @param ref_anc A \code{SummarizedExperiment} object that contains the genomic coordinates
#' and beta values of the anchor sites/regions from reference cell types.  The default
#' reference is explained in the manuscript of \pkg{decemedip}.
#' @param seed The seed for random number generation in MCMC sampling.
#' @param chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param iter A positive integer specifying the number of iterations for each chain (including
#' warmup). The default is 1000.
#' @param cores A positive integer specifying the number of cores that can be used for
#' MCMC sampling. The default is 1000.
#' @param stan_input_params A named list of parameters that specifies the prior of the
#' decemedip model.
#' @param stan_control A named list of parameters to control the sampler's behavior in
#' Stan. See the details in the documentation for the control argument in \code{\link[rstan]{stan}}.
#' @param ... Other parameters that can be passed to the \code{\link[rstan]{sampling}} function.
#'
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment ncol
#' @importFrom SummarizedExperiment nrow
#' @return
#' @export
#'
#' @examples
decemedip <- function(
    sample_bam_file = NULL,
    paired_end = NULL,
    counts_cts = c(),
    counts_anc = c(),
    ref_cts = hg19.ref.cts.se,
    ref_anc = hg19.ref.anc.se,
    diagnostics = TRUE,
    seed = 2024,
    cores = 4,
    chains = 4,
    iter = 1000,
    stan_input_params = list(
      'alpha'     = rep(1, SummarizedExperiment::ncol(ref_cts)),
      's_mu'      = 1,
      's_sigma'   = 1,
      'n_knot_z'  = 0,
      'degree_z'  = 3,
      'Xi'        = cor(as.matrix(SummarizedExperiment::assays(ref_cts)[[1]])),
      's_theta'   = 1,
      's_tau'     = 1
    ),
    stan_control = NULL,
    ...
) {

  ## Checks on reference
  stopifnot(class(ref_cts) == "RangedSummarizedExperiment")
  stopifnot(class(ref_anc) == "RangedSummarizedExperiment")
  stopifnot(SummarizedExperiment::ncol(ref_cts) == SummarizedExperiment::ncol(ref_anc))

  ## Checks on stan parameters
  stopifnot(c('alpha','s_mu','s_sigma','n_knot_z','degree_z','Xi','s_theta','s_tau')
            %in% names(stan_input_params))
  stopifnot(length(stan_input_params$alpha) == SummarizedExperiment::ncol(ref_cts))

  ## Checks on exclusivity of counts and bam files
  if (!is.null(sample_bam_file) | !is.null(paired_end)) {
    if (is.null(sample_bam_file)) stop('Please input `sample_bam_file`!')
    if (is.null(paired_end))      stop('Please input `paired_end`!')
    if (length(counts_cts) != 0)  stop('Invalid: Both counts and bam files are received as input!')
    if (length(counts_anc) != 0)  stop('Invalid: Both counts and bam files are received as input!')
  }
  if (length(counts_cts) > 0 | length(counts_anc) > 0) {
    # check whether dim aligns with reference
    stopifnot(length(counts_cts) == nrow(ref_cts))
    stopifnot(length(counts_anc) == nrow(ref_anc))
    # check whether values are integers
    stopifnot(all(counts_cts == floor(counts_cts)) && is.numeric(counts_cts))
    stopifnot(all(counts_anc == floor(counts_anc)) && is.numeric(counts_anc))
    # check on exclusivity of counts and bam files
    if (!is.null(sample_bam_file)) stop('Invalid: Both counts and bam files are received as input!')
    if (!is.null(paired_end))      stop('Invalid: Both counts and bam files are received as input!')
  }

  ## Load read counts of reference sites if bam file is provided
  if (!is.null(sample_bam_file) & !is.null(paired_end)) {
    sample_cts.se <- getRoiReadCount(sample_bam_file, sample_names = NULL, sample_paired = paired_end, roi = granges(ref_cts))
    sample_anc.se <- getRoiReadCount(sample_bam_file, sample_names = NULL, sample_paired = paired_end, roi = granges(ref_anc))
    counts_cts <- SummarizedExperiment::assays(sample_cts.se)[[1]][,1] |> unlist()
    counts_anc <- SummarizedExperiment::assays(sample_anc.se)[[1]][,1] |> unlist()
  }

  ## Check for low coverage or unusual high counts
  if (sum(counts_cts) < 0.25 * length(counts_cts))
    warning('Unusual low coverage for cell type-specific reference sites!')

  ## Prepare model variables
  y <- c(counts_cts, counts_anc)
  X <- rbind(SummarizedExperiment::assays(ref_cts)[[1]], SummarizedExperiment::assays(ref_anc)[[1]])
  z <- c(SummarizedExperiment::rowData(ref_cts)$n_cpgs_100bp, SummarizedExperiment::rowData(ref_anc)$n_cpgs_100bp) |> log1p()

  # add head and tail into knots for computation in the stan script
  if (stan_input_params$n_knot_z == 0) probs <- NULL else probs <- seq(0, 1, length.out = stan_input_params$n_knot_z)
  stan_input_params$knots_z <- quantile(z, probs = c(0, probs, 1))
  stan_input_params$n_knot_z <- length(stan_input_params$knots_z)

  ## Prepare input for stan model
  data_list <- c(list('N' = nrow(X), 'K' = ncol(X), 'y' = y, 'X' = X, 'z' = z), stan_input_params)

  ## Prepare model fitting details
  if (diagnostics) {
    model <- stanmodels$decemedip1
  } else {
    model <- stanmodels$decemedip0
  }

  if (is.null(stan_control)) stan_control <- list(adapt_delta = 0.95, max_treedepth = 15)

  ## Run MCMC
  posterior <- rstan::sampling(model, data = data_list,
                               iter = iter, chains = chains, cores = cores, seed = seed,
                               control = stan_control, ...)

  ## Return model
  return(posterior)
}
