#now we plot the results of our linear regression or linear mixed effects models
#we create 3 plots:
#1, t-statistics distribution per method per diagnostic group
#2, beta estimates distribution per method per diagnostic group
#3, barplots showing number of significant regions after FDR correction, per method per diagnostic group, and divided by effect direction

library(tidyverse)
library(ggpubr)

cols <- c('DBM\n(indirect)'='#EEBCBF',
          'DBM\n(direct)'='#A37794', 
          'Cortical\nThickness'='#9AAFBD',
          'VBM\n(indirect)'='#E2BF71', 
          'VBM\n(direct)'='#DF8673', 
          'Volume'='#BBC3AB')

#prepare input data
#(change factor level and naming for groups and methods)
prepare_data <- function(data) {
  data$EffectType <- factor(data$EffectType, levels = c('BV','SV','PNFA')) |>
    recode(BV='bvFTD', SV='svPPA', PNFA='nfvPPA')
  
  data$method <- factor(data$method, levels = c('dirDBM','DBM','VBMdir','VBMindir','CT','Vol')) |>
    recode(DBM='DBM\n(indirect)',
           dirDBM='DBM\n(direct)',
           VBMdir='VBM\n(direct)',
           VBMindir='VBM\n(indirect)',
           CT='Cortical\nThickness',
           Vol='Volume')
  
  return(data)
}

#generate plots
make_plots <- function(data, label_suffix) {
  
  ## ---- t-statistic distribution ----
  tstat <- ggplot(data, aes(y=tStat, x=method)) +
    facet_wrap(~EffectType, nrow=1) +
    geom_hline(yintercept=0, color='#3b3b3b', linetype='dashed') +
    geom_violin(aes(col=method)) +
    geom_boxplot(aes(fill=method), alpha=.6, width=.2, outliers=FALSE) +
    geom_jitter(aes(col=method), alpha=.3, width=.3) +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    ylim(-20, 15) +
    theme_minimal() +
    labs(y='regional\nt-value', x=NULL) +
    guides(col='none', fill='none') +
    theme(plot.title = element_text(size=16),
          axis.text = element_text(size=10),
          strip.text = element_text(size=12))
  
  ggsave(paste0("/data/dadmah/metame/NIFD_aim1/plots/lme2_distr_", label_suffix, "_t.png"),
         plot=tstat, dpi=600, width=12.5, height=3, units="in")
  
  
  ## ---- beta estimate distribution ----
  a <- ggplot(data, aes(y=est, x=method)) +
    facet_wrap(~EffectType, nrow=1) +
    geom_hline(yintercept=0, color='#3b3b3b', linetype='dashed') +
    geom_violin(aes(col=method)) +
    geom_boxplot(aes(fill=method), alpha=.6, width=.2, outliers=FALSE) +
    geom_jitter(aes(col=method), alpha=.3, width=.3) +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    ylim(-0.8, 0.8) +
    theme_minimal() +
    labs(y='regional\nbeta estimate', x=NULL) +
    guides(col='none', fill='none') +
    theme(plot.title = element_text(size=16),
          axis.text = element_text(size=10),
          strip.text = element_text(size=12))
  
  
  ## ---- FDR count plot ----
  fdr_counts <- data %>%
    filter(FDR == 1) %>%
    mutate(direction = ifelse(est >= 0, 1, -1)) %>%
    group_by(method, EffectType, direction) %>%
    summarise(count = n(), .groups='drop_last') %>%
    mutate(count_signed = count * direction) %>%
    group_by(method, EffectType) %>%
    mutate(
      total_count = sum(count),  # <-- total positive + negative
      total = ifelse(method == 'Cortical\nThickness', 62, 90),
      label_tot = ifelse(direction==-1,paste0(total_count, "/\n", total),''),
      label = ifelse(direction==-1 & method == 'Cortical\nThickness',paste0(count, "*"),as.character(count)),
      perc = paste0(round(total_count/total * 100,1), '%')
      ) %>%
    ungroup()
  
  b <- ggplot(fdr_counts, aes(x=method, y=count_signed, fill=method)) +
    facet_wrap(~EffectType, nrow=1) +
    geom_col(position="dodge", alpha=0.8) +
    geom_text(aes(label=label, vjust=ifelse(count_signed >= 0, -0.3, 1.2)),
              position=position_dodge(width=0.9),
              size=4.5, color="#222222") +
    geom_text(aes(label=perc, y=-110),
              position=position_dodge(width=0.9),
              size=4, color="#222222") +
    scale_y_continuous(
      limits=c(-115,30),
      breaks=seq(-100,25,by=25),
      labels=function(x) abs(x)
    ) +
    geom_hline(yintercept=0, color="#222222", linewidth=0.4,linetype="dashed") +
    scale_fill_manual(values=cols) +
    theme_minimal() +
    labs(y="significant regions after FDR correction\nnegative estimates | positive estimates",
         x=NULL) +
    guides(col='none', fill='none') +
    theme(plot.title = element_text(size=16),
          axis.text = element_text(size=10),
          strip.text = element_text(size=12),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  # Combine and save
  combined_plot <- ggpubr::ggarrange(a, b, nrow=2)
  
  ggsave(paste0("lm1_distr_", label_suffix, "_main.png"),
         plot=combined_plot, dpi=600, width=13, height=6.5, units="in")
}

#run the function
#we did all analyses in a subset (scans that passed quality control for all methods) and a full sample
#so we run it for both results here
for (suffix in c("all", "sub")) {
  file <- paste0("LME2.2_results_", suffix, "_251028.csv")
  message("Processing: ", file)
  
  data <- read.csv(file) |> prepare_data()
  make_plots(data, suffix)
}
