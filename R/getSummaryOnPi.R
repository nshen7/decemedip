#' Extract summary statistics and diagnostics on fitted cell type proportions
#'
#' @param posterior The fitted posterior object from decemedip model.
#' @param probs A numeric vector specifying quantiles of interest. The defaults
#' is c(0.025,0.25,0.5,0.75,0.975).
#' @param digits_summary The number of significant digits to use in the summary,
#' defaulting to 5.
#' @param cell_type_names Name of the cell types in reference panel. The order
#' should align with order in the reference panel.
#' @param ... Additional arguments that get passed into \code{rstan::monitor} function.
#'
#' @return A data.frame object containg summary statistics and diagnostic statistics
#' of the fitted cell type proportions.
#' @export
#'
#' @examples
#' \dontrun{
#' output <- decemedip(sample_bam_file = 'dir/to/bam/file',  paired = TRUE)
#' smr_pi.df <- getSummaryOnPi(output$posterior)
#' }
getSummaryOnPi <- function(
    posterior,
    probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
    digits_summary = 5,
    cell_type_names = NULL,
    ...
) {
  if (is.null(cell_type_names)) {
    data(hg19.ref.cts.se)
    cell_type_names = colnames(hg19.ref.cts.se)
  }
  smr_pi.df <- rstan::monitor(rstan::extract(posterior, pars = c("pi"), permuted = FALSE),
                              probs = probs,
                              digits_summary = digits_summary,
                              print = FALSE)

  if (length(cell_type_names) != nrow(smr_pi.df))
    stop('Length of cell_type_names does not align with the number of parameters posterior object.')

  smr_pi.df <- smr_pi.df |>
    as.data.frame() |>
    mutate(cell_type = factor(cell_type_names, levels = cell_type_names)) |>
    relocate(cell_type) |>
    select(1:(7+length(probs)))

  return(smr_pi.df)
}
