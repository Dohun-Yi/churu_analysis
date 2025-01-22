# bench_time_visualize

library(ggplot2)
library(ggpubr)
library(data.table)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyverse)
library(grid)
library(zoo)
library(ggExtra)
library(pheatmap)
library(gridExtra)
library(mclust)
library(factoextra)

path_root <- "../../"
path_workd <- paste(sep='',path_root, '/analysis/1_benchmark/')

# ================== Check time ================== 
path_meta <- paste(sep='',path_workd, '/time.out.txt') # all, 400k
df_meta <- fread(path_meta, skip=1, header=F)[,1:5] %>% as.data.frame()
names(df_meta) <- c('sra', 'r_star', 'r_gatk', 'r_varscan' ,'r_churu')

# wall clock
df_meta_wc <- df_meta
df_meta_wc$GATK <- df_meta_wc$r_gatk
df_meta_wc$VarScan <- df_meta_wc$r_star + df_meta_wc$r_varscan
df_meta_wc$CHURU <- df_meta_wc$r_churu

df_melt_wc <- melt(df_meta_wc, measure.vars = c('GATK','VarScan','CHURU'))
df_melt_wc$variable <- factor(df_melt_wc$variable,
                              levels=c('CHURU','VarScan','GATK'))
p1 <- ggplot(df_melt_wc, aes(x=variable, y=value/3600)) +
  geom_violin(fill='lightgrey') +
  geom_boxplot(width=0.2) +
  scale_y_log10(labels=function(x) sprintf("%.1f", x)) +
  ggtitle('Wall clock time') +
  ylab('Runtime (hours)') +
  xlab('') +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# wall clock
df_meta_cpu <- df_meta
df_meta_cpu$GATK <- (df_meta_cpu$r_gatk - df_meta_cpu$r_star) + df_meta_cpu$r_star * 80
df_meta_cpu$VarScan <- df_meta_cpu$r_star * 80 + df_meta_cpu$r_varscan
df_meta_cpu$CHURU <- df_meta_cpu$r_churu * 32

df_melt_cpu <- melt(df_meta_cpu, measure.vars = c('GATK','VarScan','CHURU'))
df_melt_cpu$variable <- factor(df_melt_cpu$variable,
                              levels=c('CHURU','VarScan','GATK'))

p2 <- ggplot(df_melt_cpu, aes(x=variable, y=value/3600)) +
  geom_violin(fill='lightgrey') +
  geom_boxplot(width=0.2) +
  scale_y_log10(labels=function(x) sprintf("%.1f", x)) +
  ggtitle('CPU time') +
  ylab('Runtime (hours)') +
  xlab('') +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

ggarrange(p1, p2, align='hv')
ggsave(filename=paste0(path_root, '/pdfs/FigureX_bench_time.pdf'), unit="in",width=4.5, height=3.0)

# ================================================ 
