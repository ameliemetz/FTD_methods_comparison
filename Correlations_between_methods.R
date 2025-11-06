#we calculated regional atrophy measures using various different methods
#here we assess how values of derived using different methods correlate with each other
#we create 3 plots:
#1, violin/boxplots showing regional correlation distribution between all methods
#2, heatmaps showing mean correlation between methods in healthy controls and FTD patients
#3, scatterplot showing correlation between TIV estimated in PELICAN vs FreeSurfer

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(forcats)

data_lme <- read.csv("LME_data_251020.csv", stringsAsFactors=TRUE)

region_cols <- tail(names(data_lme), 90)

data_lme_clean <- data_lme %>%
  filter(rowSums(!is.na(dplyr::select(., all_of(region_cols)))) > 0) %>%
  filter(visit_mri==1) %>%
  filter(DX %in% c('CON','SV','PNFA','BV'))

wide_data <- data_lme_clean %>%
  dplyr::select(DX, LONI_ID, method, all_of(region_cols)) %>%
  pivot_wider(names_from = method, values_from = all_of(region_cols), names_sep = "_")

methods <- unique(data_lme$method)

method_pairs <- combn(methods, 2, simplify = FALSE)

################################################################################
# box-/violinplots to show distribution of region-wise correlation between atrophy estimates between methods

all_corr_results <- lapply(method_pairs, function(pair) {
  m1 <- pair[1]
  m2 <- pair[2]
  
  cors <- sapply(region_cols, function(region) {
    col1 <- paste0(region, "_", m1)
    col2 <- paste0(region, "_", m2)
    x <- wide_data[[col1]]
    y <- wide_data[[col2]]
    if (all(is.na(x)) | all(is.na(y))) {
      return(NA_real_)
    }
    cor(x, y, use = "pairwise.complete.obs")
  })
  
  data.frame(region = region_cols,
             method1 = m1,
             method2 = m2,
             correlation = cors)
}) %>% bind_rows()

all_corr_results <- all_corr_results %>%
  mutate(
    method1 = case_when(
      method1 == "DBM"     ~ "DBM (indirect)",
      method1 == "dirDBM"  ~ "DBM (direct)",
      method1 == "Vol"     ~ "Volume",
      method1 == "VBMdir"  ~ "VBM (direct)",
      method1 == "VBMindir"  ~ "VBM (indirect)",
      method1 == "CT" ~ "Cortical Thickness",
      TRUE                 ~ method1
    ),
    method2 = case_when(
      method2 == "DBM"     ~ "DBM (indirect)",
      method2 == "dirDBM"  ~ "DBM (direct)",
      method2 == "Vol"     ~ "Volume",
      method2 == "VBMdir"  ~ "VBM (direct)",
      method2 == "VBMindir"  ~ "VBM (indirect)",
      method2 == "CT" ~ "Cortical Thickness",
      TRUE                 ~ method2
    )
  ) %>%
  mutate(method_pair = paste(method1, method2, sep = " - "))

all_corr_results$method_pair <- gsub(" - ", " -\n", all_corr_results$method_pair)

cor_plot <- ggplot(all_corr_results %>%
                     mutate(method_pair = fct_reorder(method_pair, correlation, .fun = median, .desc = TRUE, .na_rm = TRUE)), 
                   aes(y = correlation, x = method_pair, fill = method_pair)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed",alpha=.5) +
  geom_boxplot(alpha = .3, width = 0.1, outliers = FALSE) +
  geom_violin(aes(col = method_pair), width = 1.2, alpha = 0) +
  geom_jitter(aes(col = method_pair), width = .3, alpha = .5) +
  guides(fill = "none", color = "none") +
  colorspace::scale_color_discrete_divergingx(palette = "Zissou 1", nmax = 15) +
  colorspace::scale_fill_discrete_divergingx(palette = "Zissou 1", nmax = 15) +
  labs(#title = "Distribution of Region-wise Correlations Between Methods", 
       y = "Correlation (Pearson's r)") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("methods_corr.png", 
       plot = cor_plot, dpi = 600, width = 12, height = 3, units = "in")

################################################################################
# heatmaps, showing mean correlation between methods in controls vs FTD

all_corr_results <- lapply(method_pairs, function(pair) {
  m1 <- pair[1]
  m2 <- pair[2]
  
  cors_df <- lapply(region_cols, function(region) {
    col1 <- paste0(region, "_", m1)
    col2 <- paste0(region, "_", m2)
    
    x <- wide_data[[col1]]
    y <- wide_data[[col2]]
    dx <- wide_data$DX   # diagnosis column
    
    if (all(is.na(x)) | all(is.na(y))) {
      return(data.frame(region = region,
                        method1 = m1,
                        method2 = m2,
                        mean_cor_CON = NA_real_,
                        mean_cor_Other = NA_real_))
    }
    
    mean_cor_CON <- cor(x[dx == "CON"], y[dx == "CON"], use = "pairwise.complete.obs")
    mean_cor_Other <- cor(x[dx != "CON"], y[dx != "CON"], use = "pairwise.complete.obs")
    
    data.frame(region = region,
               method1 = m1,
               method2 = m2,
               mean_cor_CON = mean_cor_CON,
               mean_cor_Other = mean_cor_Other)
  }) %>% bind_rows()
  
  cors_df
}) %>% bind_rows()

summary_corrs <- all_corr_results %>%
  group_by(method1, method2) %>%
  summarise(
    controls = mean(mean_cor_CON, na.rm = TRUE),
    FTD = mean(mean_cor_Other, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = c(controls, FTD),
    names_to = "DX_group",
    values_to = "mean_cor"
  ) %>%
  mutate(
    method1 = case_when(
      method1 == "DBM"     ~ "DBM\n(indirect)",
      method1 == "dirDBM"  ~ "DBM\n(direct)",
      method1 == "Vol"     ~ "Volume",
      method1 == "VBMdir"  ~ "VBM\n(direct)",
      method1 == "VBMindir"  ~ "VBM\n(indirect)",
      method1 == "CT" ~ "Cortical\nThickness",
      TRUE                 ~ method1
    ),
    method2 = case_when(
      method2 == "DBM"     ~ "DBM\n(indirect)",
      method2 == "dirDBM"  ~ "DBM\n(direct)",
      method2 == "Vol"     ~ "Volume",
      method2 == "VBMdir"  ~ "VBM\n(direct)",
      method2 == "VBMindir"  ~ "VBM\n(indirect)",
      method2 == "CT" ~ "Cortical\nThickness",
      TRUE                 ~ method2
    )
  ) 

methods <- unique(c('DBM\n(direct)',"DBM\n(indirect)","VBM\n(direct)","VBM\n(indirect)","Cortical\nThickness","Volume"))

summary_corrs <- summary_corrs %>%
  rowwise() %>%
  mutate(
    m1 = ifelse(as.numeric(factor(method1, levels = methods)) <=
                  as.numeric(factor(method2, levels = methods)),
                method1, method2),
    m2 = ifelse(as.numeric(factor(method1, levels = methods)) <=
                  as.numeric(factor(method2, levels = methods)),
                method2, method1)
  ) %>%
  ungroup() %>%
  select(method1 = m1, method2 = m2, DX_group, mean_cor) %>%
  mutate(
    method1 = factor(method1, levels = methods),
    method2 = factor(method2, levels = methods)
  ) %>%
  filter(as.numeric(method1) <= as.numeric(method2))

heatdx_plot <- ggplot(summary_corrs, aes(x = method1, y = method2, fill = mean_cor)) +
  facet_wrap(~DX_group,nrow=1) +
  geom_tile(color = "white") +  
  colorspace::scale_fill_continuous_sequential(palette = "TealGrn") +
  theme_minimal() +
  geom_text(aes(label = round(mean_cor, 2)), color = "black", size = 4) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),      
    axis.line = element_blank(),      
    axis.ticks = element_blank(),      
    axis.text = element_text(color = "black")) +
  labs(fill = "Correlation (Pearson's r)", 
       #title = "Mean Region-wise Correlations Between Methods"
       )

ggsave("methods_corr_heatmap.png", 
       plot = heatdx_plot, dpi = 600, width = 10, height = 4, units = "in")

################################################################################
# TIV PELICAN vs TIV FreeSurfer

library(scales)

data_tiv <- data_lme[data_lme$DX != 'OTHER' & data_lme$visit_mri == 1,]

r_val <- cor(data_tiv$scale_TIV, data_tiv$TIV, use = "pairwise.complete.obs")
r_label <- paste0("r = ", round(r_val, 2))

tiv_corr <- ggplot(data_tiv,aes(x=scale_TIV,y=TIV,col=as.factor(GENDER))) +
  geom_point() +
  theme_minimal() +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed",alpha=0.5) +
  scale_x_continuous(limits = range(c(data_lme$scale_TIV, data_lme$TIV), na.rm = TRUE),labels = label_comma()) +
  scale_y_continuous(limits = range(c(data_lme$scale_TIV, data_lme$TIV), na.rm = TRUE),labels = label_comma()) +
  scale_color_manual(
    values = c('1' = '#66C5CC', '2' = "#0d585f"),   
    labels = c("Male", "Female")          
  ) +
  annotate("text", x = 1150000, y = 1950000,
           label = r_label, hjust = 0, vjust = 1, size = 5, color = "black") +
  labs(
    x = "PELICAN TIV",   
    y = "FreeSurfer TIV",    
    color = "Sex"                 
  )

ggsave("methods_tiv_corr.png", 
       plot = tiv_corr, dpi = 600, width = 5, height = 3, units = "in")

#############################################################################
#making a summary table

summary_table <- all_corr_results %>%
  group_by(method_pair) %>%
  summarise(
    mean_corr = mean(correlation, na.rm = TRUE),
    sd_corr = sd(correlation,na.rm=TRUE),
    range_corr = paste0(
      round(min(correlation, na.rm = TRUE), 3), " – ",
      round(max(correlation, na.rm = TRUE), 3)
    ),
    .groups = "drop"
  )
