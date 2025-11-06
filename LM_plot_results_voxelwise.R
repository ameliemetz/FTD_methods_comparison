#now we plot the voxel-/vertex-wise results of our linear regression or linear mixed effects models
#we create 3 plots:
#1, t-statistics distribution per method per diagnostic group
#2, beta estimates distribution per method per diagnostic group
#3, barplots showing number of significant regions after FDR correction, per method per diagnostic group, and divided by effect direction

library(dplyr)
library(ggplot2)
library(tidyr)

#load results files and combine PELICAN and FreeSurfer results
data_vertex <- read.csv('../results_files/LMvertex_results_sub.csv')
data_voxel <- read.csv('../results_files/voxelwise_tstat_est_fdr_251023.csv')

data_vertex$method <- 'CT'

#we created different separate figures for grey matter vs white matter results and subset vs full sample analyses
#so here we pick which one we want to plot 
data <- data_voxel %>%
  filter(cohort=='all' & method != 'VBMwm') %>%
  dplyr::select(method,group,tStat,FDR,est) %>%
  mutate(group = recode(group,
                        bvFTD = "BV",
                        svPPA = "SV",
                        nfvPPA = "PNFA"))

data <- bind_rows(
  data,
  data_vertex %>%
    dplyr::select(method, group, tStat, FDR, est) )

cols <- c('DBM\n(indirect)'='#EEBCBF','DBM\n(direct)'='#A37794', 'Cortical\nThickness'='#9AAFBD','VBM\n(indirect)'='#E2BF71', 'VBM\n(direct)'='#DF8673', 'Volume'='#BBC3AB') 

data$group <- factor(data$group, levels = c('BV','SV','PNFA'))

data$method <- factor(data$method, levels = c('dirDBM','DBM','VBM','VBMindir','CT'))

data$method <- recode(data$method, 
                      DBM='DBM\n(indirect)',
                      dirDBM='DBM\n(direct)',
                      VBM='VBM\n(direct)',
                      VBMindir='VBM\n(indirect)',
                      CT='Cortical\nThickness')

#plotting a million voxels is just too much, let's sample
data_t <- data %>%
  group_by(group, method) %>%
  slice_sample(n = 300, replace = FALSE) %>%  
  ungroup()

## ---- t-statistic distribution ----
tstat <- ggplot(data_t, aes(y=tStat, x=method)) +
  facet_wrap(~group, nrow=1) +
  geom_hline(yintercept=0, color='#3b3b3b', linetype='dashed') +
  geom_violin(aes(col=method)) +
  geom_boxplot(aes(fill=method), alpha=.6, width=.2, outliers=FALSE) +
  geom_jitter(aes(col=method), alpha=.3, width=.3) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  ylim(-20, 15) +
  theme_minimal() +
  labs(y='voxel-/vertex-wise\nt-value', x=NULL) +
  guides(col='none', fill='none') +
  theme(plot.title = element_text(size=16),
        axis.text = element_text(size=10),
        strip.text = element_text(size=12))

ggsave("lm1_distr_voxel_all_t.png",
       plot=tstat, dpi=600, width=12.5, height=3, units="in")


## ---- beta estimate distribution ----
a <- ggplot(data_t, aes(y=est, x=method)) +
  facet_wrap(~group, nrow=1) +
  geom_hline(yintercept=0, color='#3b3b3b', linetype='dashed') +
  geom_violin(aes(col=method)) +
  geom_boxplot(aes(fill=method), alpha=.6, width=.2, outliers=FALSE) +
  geom_jitter(aes(col=method), alpha=.3, width=.3) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  ylim(-5, 1.2) +
  theme_minimal() +
  labs(y='voxel-/vertex-wise\nbeta estimate', x=NULL) +
  guides(col='none', fill='none') +
  theme(plot.title = element_text(size=16),
        axis.text = element_text(size=10),
        strip.text = element_text(size=12))


## ---- FDR count plot ----
fdr_counts <- data %>%
  filter(method == 'Cortical\nThickness' & FDR == 1 | method != 'Cortical\nThickness' & FDR < 0.05) %>%
  mutate(direction = ifelse(est >= 0, 1, -1)) %>%
  group_by(method, group, direction) %>%
  summarise(count = n(), .groups='drop_last') %>%
  mutate(count_signed = count * direction) %>%
  group_by(method, group) %>%
  mutate(
    total_count = sum(count),  # <-- total positive + negative
    total = ifelse(method == 'Cortical\nThickness', 327684, 1082282),
    label_tot = ifelse(direction==-1,paste0(total_count, "/\n", total),''),
    label = ifelse(direction==-1 & method == 'Cortical\nThickness',paste0(count, "*"),as.character(count)),
    perc = paste0(round(total_count/total * 100,1), '%')
  ) %>%
  ungroup()

b <- ggplot(fdr_counts, aes(x=method, y=count_signed, fill=method)) +
  facet_wrap(~group, nrow=1) +
  geom_col(position="dodge", alpha=0.8) +
  geom_text(aes(label=label, vjust=ifelse(count_signed >= 0, -0.3, 1.2)),
            position=position_dodge(width=0.9),
            size=4.5, color="#222222") +
  geom_text(aes(label=perc, y=-900000),
            position=position_dodge(width=0.9),
            size=4, color="#222222") +
  scale_y_continuous(
    limits=c(-900000,50000),
    breaks=seq(-800000,50000,by=100000),
    labels=function(x) abs(x)
  ) +
  geom_hline(yintercept=0, color="#222222", linewidth=0.4,linetype="dashed") +
  scale_fill_manual(values=cols) +
  theme_minimal() +
  labs(y="significant voxels/vertices after FDR correction\nnegative estimates | positive estimates",
       x=NULL) +
  guides(col='none', fill='none') +
  theme(plot.title = element_text(size=16),
        axis.text = element_text(size=10),
        strip.text = element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# Combine and save
combined_plot <- ggpubr::ggarrange(a, b, nrow=2)

ggsave("lm1_distr_voxel_all_main.png",
       plot=combined_plot, dpi=600, width=13, height=6.5, units="in")

