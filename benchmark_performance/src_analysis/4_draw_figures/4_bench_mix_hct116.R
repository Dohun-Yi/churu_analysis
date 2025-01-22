# bench_mix

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
path_workd <- paste(sep='',path_root, '/analysis/3_mix_simulation/')

# ================== Check mix ==================
path_data <- paste(sep='',path_workd, '/final_result_churu.txt')
path_data_2 <- paste(sep='',path_workd, '/final_result_celid.txt')
df <- fread(path_data, header=T) %>% as.data.frame()
df2 <- fread(path_data_2, header=T) %>% as.data.frame()
dfc <- cbind (df, df2[,2:3])
dfc$mix <- str_replace_all(dfc$mix,'_',':')
df_melt <- melt(dfc, id.vars=c('mix'))
df_melt$value <- as.numeric(str_replace_all(df_melt$value,'%$',''))
df_melt$mix <- factor(df_melt$mix, levels=dfc$mix)
df_melt$legend <- ''
df_melt[df_melt$variable=='t1_d',]$legend <- 'THP-1, -expr'
df_melt[df_melt$variable=='hc_d',]$legend <- 'HCT116, -expr'
df_melt[df_melt$variable=='t1_r',]$legend <- 'THP-1, +expr'
df_melt[df_melt$variable=='hc_r',]$legend <- 'HCT116, +expr'
df_melt[df_melt$variable=='t1_c',]$legend <- 'THP-1, CeL-ID'
df_melt[df_melt$variable=='hc_c',]$legend <- 'HCT116, CeL-ID'

df_guide <- dfc
df_guide$mix <- factor(dfc$mix, levels=dfc$mix)
df_guide$t1 <- 0:10 * 10
df_guide$hc <- 10:0 * 10
my_palette <- c("t1_d" = "#1f77b4", "hc_d" = "#aec7e8",
                "t1_r" = "#d62728", "hc_r" = "#ff9896",
                "t1_c" = "black"  , "hc_c" = "grey")
my_palette <- c("d" = "#1f77b4", "r" = "#d62728", "c" = "black")

msqrt_t1_d <- floor(mean((as.numeric(str_remove_all(df_guide$t1_d, '%')) - df_guide$t1)**2, na.rm=T)*10)/10
msqrt_t1_r <- floor(mean((as.numeric(str_remove_all(df_guide$t1_r, '%')) - df_guide$t1)**2, na.rm=T)*10)/10
msqrt_t1_c <- floor(mean((as.numeric(str_remove_all(df_guide$t1_c, '%')) - df_guide$t1)**2, na.rm=T)*10)/10
msqrt_hc_d <- floor(mean((as.numeric(str_remove_all(df_guide$hc_d, '%')) - df_guide$hc)**2, na.rm=T)*10)/10
msqrt_hc_r <- floor(mean((as.numeric(str_remove_all(df_guide$hc_r, '%')) - df_guide$hc)**2, na.rm=T)*10)/10
msqrt_hc_c <- floor(mean((as.numeric(str_remove_all(df_guide$hc_c, '%')) - df_guide$hc)**2, na.rm=T)*10)/10


format(mean((as.numeric(str_remove_all(df_guide$t1_r, '%')) - df_guide$t1)**2, na.rm=T), nsmall = 5)

ggplot(df_melt, aes(x=mix)) +
  geom_line(aes(y=value, group=variable, color=legend)) +
  geom_point(aes(y=value, group=variable, color=legend, shape=legend), size=2) +
  geom_line(data=df_guide, aes(x=mix, y=t1, group=1), linetype=2, alpha=0.3) +
  geom_line(data=df_guide, aes(x=mix, y=hc, group=1), linetype=2, alpha=0.3) +
  scale_color_manual(values = c("THP-1, -expr" = "#1f77b4", "HCT116, -expr" = "#aec7e8", 
                                "THP-1, +expr" = "#d62728", "HCT116, +expr" = "#ff9896",
                                "THP-1, CeL-ID" = "black"  , "HCT116, CeL-ID" = "grey"), name='') +
  scale_shape_manual(values = c("THP-1, -expr" = 19, "HCT116, -expr" = 1, 
                                "THP-1, +expr" = 19, "HCT116, +expr" = 1,
                                "THP-1, CeL-ID" = 19, "HCT116, CeL-ID" = 1), name='') +
  annotate("text", "3:7", hjust = 0, 15, color="#1f77b4", label=paste("MSE-e=",msqrt_t1_d)) +
  annotate("text", "3:7", hjust = 0, 10, color="#d62728", label=paste("MSE+e=",msqrt_t1_r)) +
  annotate("text", "3:7", hjust = 0, 05, color="black", label=paste("MSE+e=",msqrt_t1_c)) +
  annotate("text", "3:7", hjust = 0, 100, color="#aec7e8", label=paste("MSE-e=",msqrt_hc_d)) +
  annotate("text", "3:7", hjust = 0, 095, color="#ff9896", label=paste("MSE+e=",msqrt_hc_r)) +
  annotate("text", "3:7", hjust = 0, 090, color="grey", label=paste("MSE+e=",msqrt_hc_c)) +
  xlab('Designed fraction (THP-1:HCT116)') +
  ylab('Estimated fraction') +
  # ggtitle('Modeling expression improve estimation of contamination') +
  theme_pubr() +
  theme(legend.position = 'right')

ggsave(filename=paste0(path_root, '/pdfs/FigureX_mixture.pdf'), unit="in",width=6.0, height=4.0)


# =============================================== 
