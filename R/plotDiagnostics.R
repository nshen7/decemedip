#' Diagnostics for model fitting in \code{\link{decemedip}}
#'
#' @param decemedip_output The output from \code{\link{decemedip}} function.
#' @param plot_type A string value, either 'y_fit' or 'model_fit'.
#' @param model_fit_n_samples Number of randomly selected posterior samples for
#' plotting the diagnostic plots of stan fit. For \code{plot_type = 'model_fit'}
#' only.
#' @param model_fit_label_size Label size in the plot grid. For \code{plot_type = 'model_fit'}
#' only. See the argument \code{label_size} in \code{\link[cowplot]{plot_grid}} for details.
#' @param model_fit_align Specifies how graphs in the grid should be aligned. See
#' the argument \code{align} in \code{\link[cowplot]{plot_grid}} for details.
#' @param ... Additional arguments to be fed into \code{\link[cowplot]{plot_grid}}
#' in the case of \code{plot_type = 'model_fit'}.
#'
#'
#' @return
#' @export
#'
#' @examples
plotDiagnostics <- function(
    decemedip_output,
    plot_type,
    model_fit_n_samples = 100,
    model_fit_label_size = 12,
    model_fit_align = 'hv',
    ...
) {

  data_list <- decemedip_output$data_list
  posterior <- decemedip_output$posterior

  if (! 'y_sim' %in% posterior@model_pars)
    stop('In order to get diagnostics, please run decemedip with `diagnostics = TRUE`!')

  if (plot_type != 'y_fit' & plot_type != 'model_fit')
    stop("plot_type has to be 'y_fit' or 'model_fit'!")

  y_sim <- as.matrix(posterior, pars = "y_sim")
  smr_pi.df <- rstan::monitor(
    rstan::extract(posterior, pars=c("pi"), permuted = FALSE),
    digits_summary = 5,
    print = FALSE
  )

  if (plot_type == 'y_fit') {
    plot.df <- data.frame(
      y = data_list$y,
      x = as.matrix(data_list$X) %*% matrix(smr_pi.df$mean, ncol = 1),
      z = data_list$z,
      y_pred = colMeans(y_sim),
      y_pred_2.5 = colQuantiles(y_sim, probs = 0.025),
      y_pred_97.5 = colQuantiles(y_sim, probs = 0.975)
    )

    p <- plot.df |>
      ggplot2::ggplot(ggplot2::aes(x = x)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = y_pred_2.5, ymax = y_pred_97.5), fill = 'lightgrey') +
      ggplot2::geom_linerange(aes(ymin = y, ymax = y_pred), size = 0.5, color = 'darkgrey') +
      ggplot2::geom_point(ggplot2::aes(y = y), size = 0.5, color = 'orange2') +
      ggplot2::geom_point(ggplot2::aes(y = y_pred), size = 0.5) +
      ggplot2::facet_wrap(~ exp(z), scales = 'free_y') +
      ggplot2::theme_classic()
    return(p)
  }

  if (plot_type == 'model_fit') {

    idx <- sample(nrow(y_sim), model_fit_n_samples)
    y <- data_list$y

    xmax <- quantile(y, 0.95)
    p1 <- bayesplot::ppc_dens_overlay(y, y_sim[idx, ]) +
      xlim(0, xmax) +
      ggplot2::ggtitle('Empirical density')
    p2 <- bayesplot::ppc_stat_2d(y, y_sim[idx, ], stat = c("mean", "sd")) +
      ggplot2::ggtitle('Posterior statistics (mean & sd)')
    p3 <- bayesplot::ppc_stat(y, y_sim[idx, ], stat = "max") +
      ggplot2::ggtitle('Posterior statistics (max)')
    prop_zero <- function(x) mean(x == 0)
    p4 <- bayesplot::ppc_stat(y, y_sim[idx, ], stat = "prop_zero") +
      ggplot2::ggtitle('Posterior statistics (zero proportion)')

    p <- cowplot::plot_grid(p1, p2, p3, p4, align = model_fit_align, label_size = model_fit_label_size, ...)
    return(p)
  }


}
