plotDiagnostics <- function(
    fit,
    which = c('y_fit', 'model_fit'),
    model_fit_n_samples = 100,
    model_fit_label_size = 12,
    ...
) {
  
  if (!'y_sim' %in% fit@model_pars) 
    stop('In order to get diagnostics, please run decemedip with `diagnostics = TRUE`!')
  
  y_sim <- as.matrix(fit, pars = "y_sim")
  
  if (which == 'y_fit') {
    plot.df <- data.frame(
      y = data_list$y,
      x = as.matrix(data_list$X) %*% matrix(smr_pi.df$mean, ncol = 1),
      z = data_list$z,
      y_pred = colMeans(y_sim),
      y_pred_2.5 = colQuantiles(y_sim, probs = 0.025),
      y_pred_97.5 = colQuantiles(y_sim, probs = 0.975)
    )
    
    p <- plot.df |>
      ggplot(aes(x = x)) +
      geom_ribbon(aes(ymin = y_pred_2.5, ymax = y_pred_97.5), fill = 'lightgrey') +
      geom_point(aes(y = y), size = 0.5, color = 'orange2') +
      geom_point(aes(y = y_pred), size = 0.5) +
      facet_wrap(~ exp(z), scales = 'free_y') +
      theme_classic()
  }
  
  if (which == 'model_fit') {
    idx <- sample(nrow(y_sim), model_fit_n_samples)
    
    xmax <- quantile(y, 0.95)
    p1 <- ppc_dens_overlay(y, y_sim[idx, ]) + xlim(0, xmax)
    p2 <- ppc_stat_2d(y, y_sim[idx, ], stat = c("mean", "sd"))
    p3 <- ppc_stat(y, y_sim[idx, ], stat = "max")
    prop_zero <- function(x) mean(x == 0)
    p4 <- ppc_stat(y, y_sim[idx, ], stat = "prop_zero")
    
    p <- cowplot::plot_grid(p1, p2, p3, p4, label_size = model_fit_label_size, ...)
  }
  return(p)
}
