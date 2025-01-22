# cell_heterogeneity

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
library(pROC)

path_root <- "../../"
path_workd <- paste(sep='',path_root, 'analysis/2_STR_vs_SNP')

# ================== Compare SNP and STR ===================
# STR scores
path_data_str <- paste(sep='',path_workd, '/STR_score.txt')
df_str <- fread(path_data_str, header=F) %>% as.data.frame()
names(df_str) <- c('query', 'QID', 'ref', 'name_q', 'STRsource', 'ID', 'name', 'taxa', 'Score', 'problematic', 'Scores')
df_str$key <- paste(sep='_', df_str$query, df_str$ref)
df_str$worstScore <- sapply(strsplit(df_str$Scores, ";"), function(x) min(as.numeric(x)))
df_str$bestScore <- sapply(strsplit(df_str$Scores, ";"), function(x) max(as.numeric(x)))
df_str <- unique(df_str)

cells_with_STR <- unique(df_str$query)
cells_problematic <- unique(subset(df_str, problematic==TRUE)$ref)
cells_of_interest <- setdiff(cells_with_STR, cells_problematic)

df_str <- subset(df_str, query %in% cells_of_interest & ref %in% cells_of_interest)

# Check misinforms
dup_table <- table(unique(df_str[,c('query','QID')])$query)
dup_table[dup_table>1] # KMH2, TC71, TT
# KMH2 (CVCL_1330 is correct; CVCL_S641 is not)
# TC71 (CVCL_2213 is correct; CVCL_S882 is not)
# TT (CVCL_1774 is correct; CVCL_3174 is not)
df_str <- subset(df_str, ! QID %in% c('CVCL_S641','CVCL_S882','CVCL_3174'))

nrow(df_str) # 9,090
length(cells_with_STR) # 928
length(cells_of_interest) # 866 (problematic: 62)


# SNP scores
path_data_snp <- paste(sep='',path_workd, '/matrix.txt.gz')
df_snp <- fread(path_data_snp, header=T) %>% as.data.frame()
df_snp_melt <- melt(df_snp)
names(df_snp_melt) <- c('cell1', 'cell2', 'prob')
df_snp_melt$cell1_strip <- str_remove_all(df_snp_melt$cell1, ':.*')
df_snp_melt$cell2_strip <- str_remove_all(df_snp_melt$cell2, ':.*')
df_snp_melt$key <- paste(sep='_', df_snp_melt$cell1_strip, df_snp_melt$cell2_strip)
df_snp_melt$lineage <- ifelse(str_remove_all(df_snp_melt$cell1, '.*:')==str_remove_all(df_snp_melt$cell2, '.*:'),'yes','no')
df_snp_melt$equal <- ifelse(df_snp_melt$cell1==df_snp_melt$cell2,'yes','no')
df_snp_melt$Probability <- df_snp_melt$prob * 100

length(unique(df_snp_melt$cell1_strip)) # initially 966

# remove cells without STR or problematic
df_snp_melt <- subset(df_snp_melt, cell1_strip %in% cells_of_interest & cell2_strip %in% cells_of_interest)
nrow(df_snp_melt) # 749,956
length(unique(df_snp_melt$cell1_strip)) # 866
length(unique(df_snp_melt$cell2_strip)) # 866
length(df_snp_melt$key) # 749,956
length(unique(df_snp_melt$key)) # 749,956


# merge
# df_merge <- merge(df_str, df_snp_melt, by='key')
df_merge <- merge(df_str, df_snp_melt, by='key', all.y=T)
df_melt <- melt(df_merge, measure.vars=c('Probability', 'worstScore', 'bestScore'))
df_melt$type <- 'Other'
df_melt[df_melt$lineage=='yes',]$type <- 'Lineage'
df_melt[df_melt$equal=='yes',]$type <- 'Itself'
# df_melt <- subset(df_melt, problematic == 'FALSE') # deprecated
# df_melt <- subset(df_melt, variable == 'Probability' | problematic == 'FALSE') # deprecated
df_melt$metric <- 'CHURU'
df_melt[df_melt$variable=='worstScore',]$metric <- 'STR-worst'
df_melt[df_melt$variable=='bestScore',]$metric <- 'STR-best'
df_melt$metric <- factor(df_melt$metric, levels=c('CHURU', 'STR-best', 'STR-worst'))

df_summary <- df_melt %>%
  group_by(type, metric) %>%
  summarize(n=sum(value>-1, na.rm=T),
            frac=sum(value>90, na.rm=T)/n,
            positive=sum(value>90, na.rm=T),
            negative=sum(value<90, na.rm=T))
df_summary
# 866*866 - 749068 - 22 - 866 = 0


subset(df_merge, Probability > 90 & lineage == 'no' & equal == 'no') # FP
subset(df_merge, Probability <= 90 & (lineage == 'yes' | equal == 'yes')) # FP



df_melt$value_round <- round(df_melt$value, digits=10)

# Figure X. SNP and STR score on CCLE dataset. STR similarity is drawn by CLASTR, Tanabe.
ggplot(df_melt, aes(x=metric, fill=type)) +
  geom_hline(yintercept = 90, linetype=2, alpha=0.5, color='red') +
  geom_boxplot(aes(y=value_round), color='black', alpha=0.7) +
  # geom_text(data=df_summary, aes(y=15, label=sprintf('%.2f%%',frac*100), group=type, color=type), 
  #           position=position_dodge(width = .9),angle=-70, show.legend = F) +
  # geom_point(aes(y=value, color=type), alpha=0.15, size=0.2,
  #            position=position_jitterdodge(jitter.width=0.3)) +
  scale_fill_manual(values=c('Itself'='#F05941', 'Lineage'='#872341', 'Other'='#22092C'), name='') +
  scale_color_manual(values=c('Itself'='#F05941', 'Lineage'='#872341', 'Other'='#22092C'), name='') +
  xlab('') + ylab('Match score') +
  theme_pubr() +
  theme(legend.position = 'right')

ggsave(filename=paste0(path_root, '/pdfs/FigureX_CHURUvsSTR.pdf'), unit="in",width=4.5, height=2.5)
# =====================================================
