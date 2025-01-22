# nonccle2ccle_vaf_compare
# Responsible for...
#   nonccle2ccle verification figures (main and sup)

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
library(scales)
library(lsa)
library(parallel)
library(fastglm)
library(glmnet)
library(e1071)
library(naivebayes)
library(class)
library(MASS)


path_root <- "../../"
path_var <- paste0(path_root, '/analysis/7_snp_profile_of_nonccle2ccle/parsed_variants.txt') # latest
path_sam1 <- paste0(path_root, '/analysis/7_snp_profile_of_nonccle2ccle/nonccle_records.full.txt') # latest
path_sam2 <- paste0(path_root, '/analysis/7_snp_profile_of_nonccle2ccle/related_records.full.txt') # latest
path_group <- paste0(path_root, '/analysis/7_snp_profile_of_nonccle2ccle/result_group.txt') # latest

df_var <- fread(path_var, header=F) %>% as.data.frame()
df_sam1 <- fread(path_sam1, header=F) %>% as.data.frame()
df_sam2 <- fread(path_sam2, header=F) %>% as.data.frame()
df_sam <- rbind(df_sam1, df_sam2)
df_group <- fread(path_group, header=F) %>% as.data.frame()
names(df_var) <- c('sample', 'variant', 'vaf')
names(df_sam) <- c('sample', 'call', 'claim', 'match1', 'match2', 'srr', 'proj', 'ror', 'institute')
names(df_group) <- c('pattern','sample','group','variant','vaf')
rownames(df_sam) <- df_sam$sample
df_group$group <- as.character(df_group$group)

# Get patterns of interest
patterns <- unique(df_group$pattern)
count_group <- unique(df_group[,c('pattern','sample','group')]) %>% 
  dplyr::count(pattern,group) %>%
  dcast(pattern ~ group)

names(count_group) <- c('pattern','group1','group2','group3')
count_group_target <- count_group %>%
  subset(group2 >= 2 & group3 >= 2)
count_group_skip <- count_group %>%
  subset(! pattern %in% count_group_target$pattern)

patterns_target <- count_group_target$pattern

# get Jaccard index
process_pattern_jaccard <- function(pattern_, subsample=100) {
  set.seed(1)
  df_subgroup <- subset(df_group, pattern == pattern_)
  df_info <- unique(df_subgroup[,c('sample','group')])
  df_info <- df_info %>%
    group_by(group) %>%
    slice_sample(n = subsample) %>%
    ungroup()
  
  df_subgroup <- subset(df_subgroup, sample %in% df_info$sample)
  
  variants2 <- subset(df_subgroup, group=='2') %>% dplyr::count(variant)
  variants3 <- subset(df_subgroup, group=='3') %>% dplyr::count(variant)
  
  uniqvar2 <- setdiff(subset(variants2, n>0)$variant, subset(variants3, n>0)$variant)
  uniqvar3 <- setdiff(subset(variants3, n>0)$variant, subset(variants2, n>0)$variant)
  
  variants2 <- subset(variants2, variant %in% uniqvar2)
  variants3 <- subset(variants3, variant %in% uniqvar3)
  
  variants2 <- subset(variants2, n >= 0.1 * length(unique(subset(df_subgroup, group=='2')$sample)))
  variants3 <- subset(variants3, n >= 0.1 * length(unique(subset(df_subgroup, group=='3')$sample)))
  
  if (nrow(variants2) < 3 | nrow(variants3) < 3 ){
    return(data.frame())
  }
  df_group1 <- subset(df_subgroup, group=='1')
  
  df_pred <- data.frame()
  
  for (sample_ in unique(df_group1$sample)) {
    overlap2 <- sum(subset(df_group1,sample==sample_)$variant %in% variants2$variant) / nrow(variants2)
    overlap3 <- sum(subset(df_group1,sample==sample_)$variant %in% variants3$variant) / nrow(variants3)
    df_pred <- rbind(df_pred,
                     data.frame(pattern=pattern_, sample=sample_, group='1', overlap2=overlap2, overlap3=overlap3))
  }
  return(df_pred)
}
save_pattern_heatmap3 <- function(pattern_, subsample=100, path_out=NA) {
  set.seed(1)
  df_subgroup <- subset(df_group, pattern == pattern_)
  df_info <- unique(df_subgroup[,c('sample','group')])
  df_info <- df_info %>%
    group_by(group) %>%
    slice_sample(n = subsample) %>%
    ungroup()
  df_info <- merge(df_info, df_sam, by='sample', all.x = T)
  df_info$Group <- c('1'='Cell A>B', '2'='Cell A', '3'='Cell B')[df_info$group]
  
  df_subgroup <- subset(df_subgroup, sample %in% df_info$sample)
  
  variants2 <- subset(df_subgroup, group=='2') %>% dplyr::count(variant)
  variants3 <- subset(df_subgroup, group=='3') %>% dplyr::count(variant)
  
  uniqvar2 <- setdiff(subset(variants2, n>0)$variant, subset(variants3, n>0)$variant)
  uniqvar3 <- setdiff(subset(variants3, n>0)$variant, subset(variants2, n>0)$variant)
  
  variants2 <- subset(variants2, variant %in% uniqvar2)
  variants3 <- subset(variants3, variant %in% uniqvar3)
  
  variants2 <- subset(variants2, n >= 0.1 * length(unique(subset(df_subgroup, group=='2')$sample)))
  variants3 <- subset(variants3, n >= 0.1 * length(unique(subset(df_subgroup, group=='3')$sample)))
  
  if (nrow(variants2) < 3 | nrow(variants3) < 3 ){
    print(paste('pass:', nrow(variants2), nrow(variants3)))
    return(data.frame())
  }
  
  variants <- union(variants2$variant, variants3$variant)
  
  df_subgroup <- subset(df_subgroup, sample %in% df_info$sample & variant %in% variants)
  
  df_vec <- df_subgroup %>%
    dcast(sample~variant, fill=0, value.var='vaf')
  
  row.names(df_vec) <- df_vec$sample
  row.names(df_info) <- df_info$sample
  
  df_vec <- df_vec[,-1]
  ann_colors = list(Group=c('Cell A>B'='#821131','Cell A'='#90A4D4','Cell B'='#E85C0D'))
  pheatmap(df_vec, annotation_row=df_info[,c('Group'), drop = F], border_color=NA,
           show_colnames=F, show_rownames=F, annotation_colors=ann_colors, main=pattern_,
           clustering_distance_rows='correlation', clustering_distance_cols='correlation',
           filename=path_out)
}

results <- mclapply(patterns_target, process_pattern_jaccard, mc.cores = 200) # takes some time here!
df_combined <- do.call(rbind, results)

# check failed cases
failed_patterns <- unique(patterns_target[! patterns_target %in% df_combined$pattern])
print(length(failed_patterns))
for (pattern__ in failed_patterns) {
  print(pattern__)
  print(subset(df_group, pattern == pattern__) %>% dplyr::count(group))
}

# Sort pattern
df_median <- df_combined %>%
  subset(group=='1') %>%
  group_by(pattern) %>%
  summarize(median=median(overlap3),
            mean=mean(overlap3))

df_median <- df_median[order(df_median$median),]
df_combined$pattern <- factor(df_combined$pattern,
                              levels=rev(df_median$pattern))
length(unique(df_combined_melt$pattern))

# Figure 1. Jaccard index for unique SNP set. subsample 100, SNPs >= 0.10 frequency. minSNP=3
df_combined_melt <- melt(df_combined, measure.vars=c('overlap2', 'overlap3'))
df_combined_count <- df_combined_melt %>%
  subset(variable=='overlap2') %>%
  dplyr::count(pattern)

p_base <- ggplot(df_combined_melt) +
  scale_fill_manual(values=c('overlap3'= '#E85C0D','overlap2'= '#90A4D4'),
                    labels=c('overlap3'= 'Call', 'overlap2'= 'Claim'), name='Cell') +
  scale_color_manual(values=c('overlap3'= '#821131','overlap2'= '#1E2A5E'),
                     labels=c('overlap3'= 'Call', 'overlap2'= 'Claim'), name='Cell') +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

p1 <- p_base + ylab('Jaccard index') + xlab('') +
  geom_boxplot(aes(y=value, x=pattern, color=variable, fill=variable), 
               outlier.size=0.5, show.legend = T) +
  theme(axis.text.x = element_blank())

p2 <- p_base + ylab('nSample') + xlab('') +
  geom_bar(data=df_combined_count, stat='identity', aes(x=pattern, y=n))

ggarrange(p1, p2, align='v', nrow=2, ncol=1, heights=c(3,3))
ggsave(filename=paste0(path_root, '/pdfs/20241111_FigureX_nonccle2ccle_jaccard_index.pdf'), unit="in",width=12.0, height=6.0)


# Sup Figures. heatmap for all patterns. Clustering by correlation
drawn <- list()
for (i in seq_along(df_median$pattern)) {
  pattern_ <- rev(df_median$pattern)[i]
  pattern_name <- str_replace_all(pattern_, '/',' ')
  path_out_ <- paste0(path_root, "/pdfs/20240926_SupFigs_heatmaps/SupplementaryFigure", i, "_", pattern_name,".pdf")
  path_out_ <- str_replace_all(path_out_, ' > ','_to_')
  print(path_out_)
  drawn <- c(drawn, pattern_)
  save_pattern_heatmap3(pattern_, path_out=path_out_)
}


# checking some numbers
nrow(df_sam1)
nrow(subset(df_sam1,V5=='contam-unknown-nonccle')) # 1,017
nrow(df_median) # analised patterns: 85
nrow(subset(df_combined,group==1)) # analised samples: 498 (unanalised: 519)
