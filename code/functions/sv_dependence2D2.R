# Recode to allow use of divergent palette
sv_dependence2D2<- function(object, x, y,
                            viridis_args = getOption("shapviz.brewer_args"),
                            jitter_width = NULL, 
                            jitter_height = NULL,
                            interactions = FALSE, 
                            add_vars = NULL, ...) {
  p <- max(length(x), length(y))
  if (p > 1L) {
    if (is.null(jitter_width)) {
      jitter_width <- replicate(p, NULL)
    }
    if (is.null(jitter_height)) {
      jitter_height <- replicate(p, NULL)
    }
    plot_list <- mapply(
      FUN = sv_dependence2D2,
      x = x,
      y = y,
      jitter_width = jitter_width,
      jitter_height = jitter_height,
      MoreArgs = list(
        object = object,
        viridis_args = viridis_args,
        interactions = interactions,
        ...
      ),
      SIMPLIFY = FALSE
    )
    return(patchwork::wrap_plots(plot_list))
  }
  
  S <- get_shap_values(object)
  X <- get_feature_values(object)
  S_inter <- get_shap_interactions(object)
  nms <- colnames(object)
  stopifnot(
    x %in% nms,
    y %in% nms,
    is.null(add_vars) || all(add_vars %in% nms)
  )
  if (interactions && is.null(S_inter)) {
    stop("No SHAP interaction values available in 'object'.")
  }
  
  # Set jitter value
  if (is.null(jitter_width)) {
    jitter_width <- 0.2 * .is_discrete(X[[x]], n_unique = 7L)
  }
  if (is.null(jitter_height)) {
    jitter_height <- 0.2 * .is_discrete(X[[y]], n_unique = 7L)
  }
  
  # Color variable
  if (!interactions) {
    s <- rowSums(S[, c(x, y, add_vars)])
  } else {
    s <- S_inter[, x, y]
    if (x != y) {
      s <- 2 * s  # Off-diagonals need to be multiplied by 2 for symmetry reasons
    }
  }
  dat <- data.frame(SHAP = s, X[, c(x, y)], check.names = FALSE)
  vir <- ggplot2::scale_color_gradient2
  if (is.null(viridis_args)) {
    viridis_args <- list()
  }
  my_plot <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[x]], y = .data[[y]], color = SHAP)) +
    ggplot2::geom_jitter(width = jitter_width, height = jitter_height, ...) +
    do.call(vir, viridis_args) +
    theme_classic() +
    ggplot2::theme(legend.box.spacing = grid::unit(0, "pt"),
                   axis.ticks = element_blank(),
                   plot.title = element_text(size = 19, family = "serif"),
                   axis.text = element_text(family = "serif", size = 16),
                   axis.title = element_text(family = "serif", size = 19),
                   axis.text.x = element_text(angle = 45, vjust = 0.7),
                   strip.text = element_text(family = "serif", size = 19),
                   legend.title = element_text(family = "serif", size = 17),
                   legend.text = element_text(family = "serif", size = 16))
  change_legend_breaks(my_plot, "colour", 
                       breaks = c(min(my_plot$plot_env$s), 0, max(my_plot$plot_env$s)),
                       labels = c("low", "0", "high")) # change the legend labels
  return(my_plot)
}
