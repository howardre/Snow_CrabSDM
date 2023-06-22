rel_inf <- function(base_brt, title){
  windows()
  effects <- tibble::as_tibble(summary.gbm(base_brt$model, plotit = FALSE))
  effects %>% arrange(desc(rel.inf)) %>%
    ggplot(aes(x = forcats::fct_reorder(.f = var,
                                        .x = rel.inf),
               y = rel.inf,
               fill = rel.inf)) +
    geom_col() +
    coord_flip() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = 'Variable',
         y = 'Relative Influence',
         title = title) +
    theme_minimal() +
    theme(legend.position = "none",       
          plot.title = element_text(size = 22, family = "serif", face = "bold"),
          axis.text = element_text(family = "serif", size = 16),
          axis.title = element_text(family = "serif", size = 20),
          strip.text = element_text(family = "serif", size = 20)) 
}