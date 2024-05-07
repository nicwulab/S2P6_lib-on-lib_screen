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

plot_bind_vs_exp <- function(df, graphname){
  print (paste('correlation for:', graphname, cor(df$avg_exp_score, df$avg_norm_binding_score)))
  print (paste('# of variants:', length(df$avg_exp_score)))
  textsize <- 8
  p <- ggplot(df,aes(x=avg_exp_score, y=avg_norm_binding_score)) +
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
    labs(x='Expression score',y='Binding score')
  ggsave(graphname, p, height=2, width=2, dpi=600)
}

plot_histogram <- function(df, xmin, xmax, xlab, graphname){
  textsize <- 8
  p <- ggplot(df,aes(x=param)) +
    geom_histogram(binwidth=0.1, fill='grey30', alpha=0.3, color='black') +
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
    xlim(xmin, xmax) +
    labs(x=xlab,y='Count')
  ggsave(graphname, p, height=2, width=2, dpi=600)
  }

df <- read_tsv('result/mut_scores.tsv') %>%
        filter(totalfreq > 0.00001)
plot_bind_vs_exp(df, 'graph/bind_vs_exp.png')

plot_histogram(mutate(df,param=avg_exp_score), 0, 1.1, 'Expression score', 'graph/hist_exp.png')
plot_histogram(mutate(df,param=avg_norm_binding_score), -1.2, 1.5, 'Binding score', 'graph/hist_bind.png')
