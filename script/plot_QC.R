#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(data.table)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(gridExtra)
library(viridis)
library(qualpalr)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_replicate_cor <- function(df, graphname){
  print (paste('correlation for:', graphname, cor(df$rep1, df$rep2)))
  print (paste('# of variants:', length(df$rep1)))
  textsize <- 7
  p <- ggplot(df,aes(x=rep1, y=rep2)) +
    geom_point(size=0.1, alpha=0.3, color='grey30', pch=16) +
    theme_cowplot(12) +
    theme(plot.title=element_text(size=textsize,face="bold",hjust=0.5),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_text(size=textsize,face="bold"),
          legend.text=element_text(size=textsize,face="bold"),
          legend.position='right') +
    labs(x='Binding score (replicate 1)',y='Binding score (replicate 2)')
  ggsave(graphname, p, height=2, width=2, dpi=600)
}

plot_freq_cutoff_vs_cor <- function(df, graphname, ylab){
  textsize <- 7
  p <- ggplot(df,aes(x=freqs, y=param)) +
    geom_point(size=1, alpha=1, color='black', pch=16) +
    theme_cowplot(12) +
    theme(plot.title=element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_text(size=textsize,face="bold"),
          legend.text=element_text(size=textsize,face="bold"),
          legend.position='right') +
    labs(x='cutoff frequency',y=ylab)
  ggsave(graphname, p, height=3, width=3, dpi=600)
}

freq_analysis <- function(df){
  freqs <- seq(0,0.0001,0.000005)
  cors  <- c()
  n <- c()
  for (freq in freqs){
    df_filter <- df %>% filter(totalfreq > freq)
    cors <- c(cors, cor(df_filter$Rep1_norm_binding_score, df_filter$Rep2_norm_binding_score))
    n <- c(n, length(df_filter$totalfreq))
    }
  return (data.table(cbind(freqs, cors, n)))
  }

df <- read_tsv('result/mut_scores.tsv')
df_freq <- freq_analysis(df)
plot_freq_cutoff_vs_cor(mutate(df_freq, param=n), 'graph/QC_cutoff_freq_vs_variant_num.png', '# of variants')
plot_freq_cutoff_vs_cor(mutate(df_freq, param=cors), 'graph/QC_cutoff_freq_vs_cor.png', 'correlation between replicates')

df <- df %>%
        filter(totalfreq > 0.00001) %>%
        mutate(rep1=Rep1_norm_binding_score, rep2=Rep2_norm_binding_score)
plot_replicate_cor(df, 'graph/QC_replicate_cor.png')
