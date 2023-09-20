library(ggplot2)
library(here)
library(dplyr)
library(reshape2)

crab_summary <- readRDS(here('data/Snow_CrabData', 'crab_pca.rds')) %>%
  dplyr::select(-geometry)

# Transform female and male data
crab_trans <- mutate(crab_summary,
                     lncount_mat_female = log(mature_female + 1),
                     lncount_imm_female = log(immature_female + 1),
                     lncount_leg_male = log(legal_male + 1),
                     lncount_sub_male = log(sublegal_male + 1),
                     pres_imm_female = ifelse(immature_female > 0, 1, 0),
                     pres_mat_female = ifelse(mature_female > 0, 1, 0),
                     pres_leg_male = ifelse(legal_male > 0, 1, 0),
                     pres_sub_male = ifelse(sublegal_male > 0, 1, 0),
                     year_f = as.factor(year),
                     log_pcod_cpue = log(pcod_cpue + 1)) %>%
  filter(!is.na(temperature),
         !is.na(julian),
         !is.na(depth),
         !is.na(ice_mean),
         year_f != 2022) 


# Create plots of average CPUE over time per stage/sex
cpue_sums <- crab_trans %>%
  group_by(year) %>%
  summarise("Mature Female" = sum(mature_female, na.rm = TRUE),
            "Immature Female" = sum(immature_female, na.rm = TRUE),
            "Legal Male" = sum(legal_male, na.rm = TRUE),
            "Sublegal Male" = sum(sublegal_male, na.rm = TRUE))

cpue_tall <- melt(cpue_sums, 
                  id.vars = "year")

ggplot(cpue_tall, 
       aes(x = year,
           y = value)) +
  geom_area(color = "darkred", 
            size = 1,
            fill = "darkred",
            alpha = 0.4) +
  geom_hline(yintercept = 0) +
  facet_grid(variable ~ .) +
  theme_minimal() +
  labs(x = "Year",
       y = "Abundance") +
  theme(axis.ticks = element_blank(),
        plot.title = element_text(size = 17,
                                  family = "serif",
                                  face = "bold"),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0, "null"),
        axis.title = element_text(family = "serif", size = 17),
        axis.text.x = element_text(angle = 0, vjust = 0.7),
        strip.text = element_text(family = "serif", size = 17))

