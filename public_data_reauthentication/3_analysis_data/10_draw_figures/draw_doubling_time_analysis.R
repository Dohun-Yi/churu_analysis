# doubling_time_analysis
# Responsible for...
#   doubling time analysis (log2 fold change and by occurrences)

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
library(ggrepel)
library(tools)

path_root <- "../../"
path_dt <- paste0(path_root, '/analysis/9_doubling_time/doubling_time_fold_change.txt') # latest
path_institute <- paste0(path_root, '/analysis/9_doubling_time/institute_selected.parsed.corrected.add_time.txt') # might be obsolete in this version, won't fix
path_table <- paste0(path_root, '/analysis/9_doubling_time/final_table.txt') # latest

df_dt <- fread(path_dt, header=F) %>% as.data.frame()
df_institute <- fread(path_institute, header=F) %>% as.data.frame()
df_table <- fread(path_table, header=T) %>% as.data.frame()

df_merge <- merge(df_dt, df_institute[,c('V1','V3','V5')], by='V1', all.x=T)
df_merge <- merge(df_merge, df_table[,c('biosample','bioproject')], by.x='V1', by.y='biosample', all.x=T)
names(df_merge) <- c('sam', 'call', 'claim', 'curation', 'call_dt', 'claim_dt', 'log2fc', 'ROR_ID', 'country', 'bioproject')

nrow(df_merge) # 203,331
nrow(subset(df_merge, curation=='contam-unknown-ccle')) # 3,314

df_merge_filt <- df_merge[df_merge$log2fc != 'N/A',]
df_merge_filt$log2fc <- as.numeric(df_merge_filt$log2fc)

nrow(df_merge_filt) # 71,248



# ============= Does the contaminants has relatively shorter doubling time? ===============
df_merge_filt <- df_merge_filt[order(df_merge_filt$log2fc),]
df_merge_filt$sam <- factor(df_merge_filt$sam,
                            levels=df_merge_filt$sam)
df_ccle_contam <- subset(df_merge_filt, curation=='contam-unknown-ccle')

df1 <- unique(df_ccle_contam[,c('call','claim','log2fc')])

nrow(df_ccle_contam) # 2,761
nrow(df1) # 415


# Figure. on 415 unknown CCLE contamination patterns, contaminant has shorter doubling time in 243 (58.6%) cases.
ggplot(df1) +
  geom_density(aes(x=log2fc), fill='#666666', alpha=0.5) +
  geom_vline(xintercept=0, linetype=2, color='black', alpha=0.5) +
  ylab('Density') +
  xlab('Log2 fold change in doubling time (Call/Claim)')+
  theme_pubr()

ggsave(filename=paste0(path_root, '/pdfs/20241024_FigureX_doubling_time_fold_change_pdf.pdf'), unit="in", width=5.0, height=3.0)



# ============= Frequent contaminant cells has shorter doubling time? - bioproject level ===============
process_dt <- function(df, cut=10) {
  patterns_count <- dplyr::count(df, curation, cell) %>% dcast(cell~curation, fill=0)
  patterns_count$percent <- patterns_count$`contam-unknown-ccle` / (patterns_count$`contam-unknown-ccle` + patterns_count$match)
  patterns_count$total <- patterns_count$`contam-unknown-ccle` + patterns_count$match
  patterns_count <- merge(patterns_count, unique(df[,c('cell','dt')]), all.x=T)
  patterns_count$bin <- "0"
  patterns_count[patterns_count$`contam-unknown-ccle` > 0,]$bin <- "1"
  patterns_count[patterns_count$`contam-unknown-ccle` > 1,]$bin <- "2"
  patterns_count[patterns_count$`contam-unknown-ccle` > 2,]$bin <- "3-"
  patterns_count$bin <- factor(patterns_count$bin, 
                               levels=c("0","1","2","3-"))
  return(subset(patterns_count, total >= cut))
}
ttest_pval <- function(df) {
  dt0 <- subset(df,bin=="0")$dt
  dt3 <- subset(df,bin=="3-")$dt
  return(t.test(dt0, dt3)$p.value)
}

df_call_all <- unique(subset(df_merge_filt, curation %in% c('contam-unknown-ccle','match'))[,c('call','curation','call_dt','bioproject')])
names(df_call_all) <- c('cell','curation','dt','bioproject')

df_call_all_count <- process_dt(df_call_all, 10)

length(unique(df_call_all$cell)) # 881
nrow(df_call_all_count) # 169

# Figure. Doubling time comparison for prevalent contaminant, n=169, for publication
set.seed(5)
ggplot(df_call_all_count) +
  geom_boxplot(aes(y=dt, x=bin, group=bin, fill=bin), outliers=F, alpha=0.5, show.legend = F) +
  geom_jitter(aes(y=dt, x=bin, group=bin, color=bin), alpha=1, show.legend = F) +
  geom_segment(aes(x=1, xend=4, y=max(dt)*1.1, yend=max(dt)*1.1), size=0.3) +
  geom_segment(aes(x=1, xend=1, y=max(dt)*1.1, yend=max(dt)*1.06), size=0.3) +
  geom_segment(aes(x=4, xend=4, y=max(dt)*1.1, yend=max(dt)*1.06), size=0.3) +
  scale_y_log10(breaks=c(1:15*20)) +
  scale_color_manual(values=c("0"="#666666", "1"="#A38370","2"="#D24D3B", "3-"="#B31312")) +
  scale_fill_manual(values=c("0"="#666666", "1"="#A38370","2"="#D24D3B", "3-"="#B31312")) +
  xlab('Misidentification occurrences') + ylab('Doubling time in hours') +
  theme_pubr() +
  annotate("text", x=-Inf, y=Inf, hjust=-0.5, vjust=1, size=4,
           label=paste0("italic(p) ==", signif(ttest_pval(df_call_all_count), digits=3)), parse=TRUE)
ggsave(filename=paste0(path_root, '/pdfs/20241024_FigureX_doubling_time_by_occurrences.pdf'), unit="in", width=4.0, height=4.0)

# Check some points
df_call_all_count[1,]
df_call_all_count %>% slice_max(dt, n=5)
df_call_all_count %>% slice_min(dt, n=5)
df_call_all_count %>% slice_max(total, n=20)
