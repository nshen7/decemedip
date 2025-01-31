#' Extract summary statistics and diagnostics on fitted cell type proportions
#'
#' @param posterior The fitted posterior object from decemedip model.
#' @param probs A numeric vector specifying quantiles of interest. The defaults
#' is c(0.025,0.25,0.5,0.75,0.975).
#' @param digits_summary The number of significant digits to use in the summary,
#' defaulting to 5.
#' @param ... Additional arguments that get passed into \code{rstan::monitor} function.
#'
#' @return
#' @export
#'
#' @examples
getSummaryOnPi <- function(
    posterior,
    probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
    digits_summary = 5,
    ...
) {
  smr_pi.df <-
    rstan::monitor(rstan::extract(posterior, pars = c("pi"), permuted = FALSE),
                   probs = probs,
                   digits_summary = digits_summary,
                   print = FALSE) |>
    as.data.frame() |>
    mutate(cell_type = factor(colnames(hg19.ref.cts.se), levels = colnames(hg19.ref.cts.se))) |>
    relocate(cell_type) |>
    select(1:(7+length(probs)))

  return(smr_pi.df)
}
