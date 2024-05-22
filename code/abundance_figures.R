library(ggplot2)
library(here)
library(dplyr)
library(reshape2)
library(readr)

crab_abun <- read_csv(here('data', 'snowcrab_abundance.csv'))

# Create plots of average CPUE over time per stage/sex
crab_tall <- melt(crab_abun, 
                  id.vars = "Year")
crab_tall$variable <- as.character(crab_tall$variable)

abun_final <- subset(crab_tall,
                     !(endsWith(variable, 'crab_tall') | endsWith(variable, '_ci')))

ggplot(abun_final, 
       aes(x = Year,
           y = value)) +
  geom_area(color = "darkcyan", 
            size = 1,
            fill = "darkcyan",
            alpha = 0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(factor(variable, levels = c("Legal Male", 
                                         "Sublegal Male", 
                                         "Mature Female", 
                                         "Immature Female")) ~ .,
             scales = "free_y") +
  theme_minimal() +
  labs(y = "Abundance (millions)") +
  theme(axis.ticks = element_blank(),
        plot.title = element_text(size = 17,
                                  family = "serif",
                                  face = "bold"),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0, "null"),
        axis.title = element_text(family = "serif", size = 17),
        axis.text.x = element_text(angle = 0, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 17))
