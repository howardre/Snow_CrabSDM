sv_dependence2 <- function(object, v, color_var = "auto", color = "#3b528b",
                                  viridis_args = getOption("shapviz.viridis_args"),
                                  jitter_width = NULL, interactions = FALSE, ...) {
  p <- length(v)
  if (p > 1L || length(color_var) > 1L) {
    if (is.null(color_var)) {
      color_var <- replicate(p, NULL)
    }
    if (is.null(jitter_width)) {
      jitter_width <- replicate(p, NULL)
    }
    plot_list <- mapply(
      FUN = sv_dependence2,
      v = v,
      color_var = color_var,
      color = color,
      jitter_width = jitter_width,
      MoreArgs = list(
        object = object,
        viridis_args = viridis_args,
        interactions = interactions,
        ...
      ),
      SIMPLIFY = FALSE
    )
    nms <- if (length(v) > 1L) v
    plot_list <- add_titles(plot_list, nms = nms)  # see sv_waterfall()
    return(patchwork::wrap_plots(plot_list))
  }
  
  S <- get_shap_values(object)
  X <- get_feature_values(object)
  S_inter <- get_shap_interactions(object)
  nms <- colnames(object)
  stopifnot(
    v %in% nms,
    is.null(color_var) || (color_var %in% c("auto", nms))
  )
  if (interactions && is.null(S_inter)) {
    stop("No SHAP interaction values available in 'object'.")
  }
  
  # Set jitter value
  if (is.null(jitter_width)) {
    jitter_width <- 0.2 * .is_discrete(X[[v]], n_unique = 7L)
  }
  
  # Set color value
  if (!is.null(color_var) && color_var == "auto" && !("auto" %in% nms)) {
    scores <- potential_interactions(object, v)
    color_var <- names(scores)[1L]  # NULL if p = 1L
  }
  if (isTRUE(interactions)) {
    if (is.null(color_var)) {
      color_var <- v
    }
    if (color_var == v) {
      y_lab <- "SHAP main effect"
    } else {
      y_lab <- "SHAP interaction"
    }
    s <- S_inter[, v, color_var]
    if (color_var != v) {
      s <- 2 * s  # Off-diagonals need to be multiplied by 2 for symmetry reasons
    }
  } else {
    y_lab <- "SHAP value"
    s <- S[, v]
  }
  dat <- data.frame(s, X[[v]])
  colnames(dat) <- c("shap", v)
  if (is.null(color_var) || color_var == v) {
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[v]], y = shap)) +
      ggplot2::geom_jitter(color = color, width = jitter_width, height = 0, ...) +
      ggplot2::ylab(y_lab)
    return(p)
  }
  dat[[color_var]] <- X[[color_var]]
  if (.is_discrete(dat[[color_var]], n_unique = 0L)) {  # only if non-numeric
    vir <- ggplot2::scale_color_viridis_d
  } else {
    vir <- ggplot2::scale_color_viridis_c
  }
  if (is.null(viridis_args)) {
    viridis_args <- list()
  }
  ggplot2::ggplot(
    dat, ggplot2::aes(x = .data[[v]], y = shap, color = .data[[color_var]])
  ) +
    ggplot2::geom_jitter(width = jitter_width, height = 0, ...) +
    ggplot2::ylab(y_lab) +
    do.call(vir, viridis_args) +
    theme_classic() +
    ggplot2::theme(legend.box.spacing = grid::unit(0, "pt"),
                   axis.ticks = element_blank(),
                   plot.title = element_text(size = 22, family = "serif", face = "bold"),
                   axis.text = element_text(family = "serif", size = 18),
                   axis.title = element_text(family = "serif", size = 22),
                   axis.text.x = element_text(angle = 0, vjust = 0.7),
                   strip.text = element_text(family = "serif", size = 21),
                   legend.title = element_text(family = "serif", size = 18),
                   legend.text = element_text(family = "serif", size = 17))
}

# Helper functions
.is_discrete <- function(z, n_unique) {
  is.factor(z) || is.character(z) || is.logical(z) || (length(unique(z)) <= n_unique)
}

# Bins z into integer valued bins, but only if discrete
.fast_bin <- function(z, n_bins) {
  if (.is_discrete(z, n_unique = n_bins)) {
    return(z)
  }
  q <- stats::quantile(z, seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE)
  findInterval(z, unique(q), rightmost.closed = TRUE)
}