rel_inf_fig <- function(inf_list){
  all_rel_inf <- Reduce(merge, 
                        lapply(inf_list,
                               function(x) 
                                 data.frame(x, rel_inf = row.names(x))))
  rownames(all_rel_inf) <- all_rel_inf[, 1]
  colnames(all_rel_inf) <- c("vars", "mature female", "immature female", "legal male", "sublegal male")
  rownames(all_rel_inf)[rownames(all_rel_inf) == "ice_mean"] <- "ice"
  rownames(all_rel_inf)[rownames(all_rel_inf) == "log_pcod_cpue"] <- "cod CPUE"
  rownames(all_rel_inf)[rownames(all_rel_inf) == "julian"] <- "day of year"
  ord_rel_inf <- all_rel_inf %>%
    arrange(factor(vars,
                   levels = c("longitude", "depth", "phi", "temperature", 
                              "latitude", "julian", "ice_mean",
                              "log_pcod_cpue", "fishery loading", "BCS")))
  mat_rel_inf <- as.matrix(ord_rel_inf[, -1])
  sca_rel_inf <- apply(mat_rel_inf,
                       MARGIN = 2,
                       FUN = function(X) 
                         (X - min(X))/diff(range(X))) # scale from 0 to 1 for each column
  heatmap_palette <- colorRampPalette(c("lightgray", "deeppink4"))
  heatmap.2(sca_rel_inf, 
            dendrogram = "none",
            Rowv = FALSE,
            Colv = FALSE,
            col = heatmap_palette(100),
            cexCol = 1.2,
            # key = TRUE,
            # density.info = "none",
            # key.title = "Relative Influence",
            # key.ytickfun = FALSE,
            # key.xlab = "",
            trace = "none",
            srtCol = 45,
            margins = c(12, 8),
            adjCol = c(NA, -0.2),
            offsetCol = 0)
}
