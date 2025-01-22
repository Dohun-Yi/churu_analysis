# ccle2ccle_vaf_compare
# Responsible for...
#   CCLE to CCLE contamination figures (hemogenous, heterogenous, and supplementary)

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
library(ggbeeswarm)
library(ggbreak)

path_root <- "../../"
path_var <- paste0(path_root, '/analysis/6_snp_profile_of_ccle2ccle/parsed_variants.txt') # latest
path_ref <- paste0(path_root, '/analysis/6_snp_profile_of_ccle2ccle/reference_snp_numbers.txt')
path_info <- paste0(path_root, '/analysis/8_final_tables/final_table.txt')
df_var_ori <- fread(path_var, header=F) %>% as.data.frame()
df_ref <- fread(path_ref, header=F) %>% as.data.frame()
df_info <- fread(path_info, header=T) %>% as.data.frame()

df_var <- df_var_ori

names(df_var) <- c('sample', 'cell', 'name', 'ach', 'vaf_obs', 'vaf_ref', 'genotype', 'frac', 'source', 'mode')
df_var$group <- paste0(df_var$cell, " - ", df_var$genotype)

names(df_ref) <- c('ach', 'nSNP')


###### Sort data ###### 
# get median VAF of call hetero
df_median <- df_var %>% 
  subset(group == 'call - 1') %>%
  group_by(sample) %>% 
  summarize(median = median(vaf_obs))

df_mean_call <- df_var %>% 
  subset(cell == 'call') %>%
  group_by(sample) %>% 
  summarize(mean = mean(vaf_obs))

# get churu estimate
df_frac <- unique(df_var[,c('sample', 'cell', 'frac')])
df_frac <- merge(df_frac, df_median, by='sample')
df_frac <- merge(df_frac, df_mean_call, by='sample')
df_frac_subset <- subset(df_frac, cell=='call')

# sort order
setorder(df_frac_subset, -frac) # by frac

df_var$order <- factor(df_var$sample, levels=df_frac_subset$sample)
df_frac$order <- factor(df_frac$sample, levels=df_frac_subset$sample)


###### Plotting - boxplot seperate for misidentification and cross-contamination ###### 
# Get samples with Claim SNP
df_n <- df_var %>% 
  subset(cell == 'claim') %>%
  dplyr::count(order)

length(unique(df_var$sample)) # 3,314 (total CCLE2CCLE samples)
sum(df_n$n >= 1) # 1,048
sum(df_n$n >= 2) # 476
sum(df_n$n >= 3) # 268

samples_cross <- subset(df_n, n>=3)$order; length(samples_cross) # 268
samples_misid <- unique(subset(df_var, ! order %in% samples_cross)$order); length(samples_misid) # 3,046


# Cross contamination (heterogenous)
df_var_cross <- subset(df_var, order %in% samples_cross)
df_var_cross_mixture <- unique(df_var_cross[,c('order','cell','frac')])
df_var_cross_stat <- df_var_cross %>%
  group_by(order, group, cell) %>%
  summarize(ymin = quantile(vaf_obs, 0.00),
            lower = quantile(vaf_obs, 0.25),
            median = median(vaf_obs),
            upper = quantile(vaf_obs, 0.75),
            ymax = quantile(vaf_obs, 1.00),
            total=n())
df_count <- df_var_cross %>%
  dplyr::count(order, cell)

mypalette <- c('call - 1'= '#E85C0D', 'call - 2'= '#821131',
               'claim - 1'= '#90A4D4', 'claim - 2'= '#1E2A5E')
mylabels <- c('call - 1'= 'Hetero', 'call - 2'= 'Homo',
                 'claim - 1'= 'Hetero', 'claim - 2'= 'Homo')

draw_it <- function(cell_type) { # call or claim
  # subset from global variables
  tmp_df_var_cross <- subset(df_var_cross, cell==cell_type)
  tmp_df_var_cross_stat <- subset(df_var_cross_stat, cell==cell_type)
  tmp_df_var_cross_mixture <- subset(df_var_cross_mixture, cell==cell_type)
  tmp_df_count <- subset(df_count, cell==cell_type)
  
  plot_base <- ggplot(tmp_df_var_cross) +
    scale_fill_manual(values=mypalette, labels=mylabels, name='genotype') +
    scale_color_manual(values=mypalette, labels=mylabels, name='genotype') +
    theme_pubr() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),)
  
  plot_VAF <- plot_base +
    geom_hline(yintercept=c(0, 0.5, 1), linetype=2, alpha=0.25) +
    geom_boxplot(data=tmp_df_var_cross,
                 outlier.shape = NA, position = position_dodge2(preserve = "single"), lwd=0.5,
                 aes(x=order, y=vaf_obs, group = interaction(order, group),
                     fill=group, color=group), alpha=1.0) +
    geom_boxplot(data=tmp_df_var_cross_stat, stat='identity',
                 outlier.shape = NA, position = position_dodge2(preserve = "single"), # This is for median
                 aes(ymin=median, ymax=median, # No wiskers
                     x=order, lower=median, middle=median, upper=median, group = interaction(order, group),
                 ), color='white', fill=NA, lwd=0.5, show.legend = F) +
    geom_line(data=tmp_df_var_cross_mixture, aes(x = order, y = frac, group=1), color='red', linetype=1) +
    ylab('Varient allele frequency') + xlab('') +
    theme(axis.title.x=element_blank())
  
  plot_nSNP <- plot_base + 
    geom_hline(yintercept=c(10,100,1000), linetype=2, alpha=0.25) +
    geom_bar(data=tmp_df_count, stat='identity', width=1,
             aes(x=order, y=n)) +
    scale_y_log10(limits=c(1,5000)) +
    ylab('nSNP') + xlab('Samples')
  
  return (ggarrange(plot_VAF, plot_nSNP, ncol=1, nrow=2, align='v', heights=c(4,2.5), common.legend = T, legend='right'))
}

# Figure. CCLE2CCLE VAF by boxplot, cross-contamination only (n=268). claim nSNP >= 3. - order (-frac). top=Call, bottom=Claim.
p1 <- draw_it('call')
p2 <- draw_it('claim')
ggarrange(p1, p2, ncol=1, nrow=2)
ggsave(filename=paste0(path_root, '/pdfs/20240910_FigureX_ccle2ccle_vaf_by_boxplot_cross.pdf'), unit="in",width=12.0, height=6.0)


###### Count by Fisher's exact test 20241106 #####
df_count_cross <- dplyr::count(df_var_cross, sample, cell, ach)
df_count_cross <- merge(df_count_cross, df_ref, by='ach', all.x=T)
df_count_table <- merge(subset(df_count_cross, cell=='call'),
                        subset(df_count_cross, cell=='claim')[,c('sample','n','nSNP')], 
                        by='sample')

names(df_count_table) <- c('sample','ach','v1','call_found','call_all','claim_found','claim_all')
df_count_table$call_notfound <- df_count_table$call_all - df_count_table$call_found
df_count_table$claim_notfound <- df_count_table$claim_all - df_count_table$claim_found

nrow(df_count_table) # 268

df_count_table$pvalue <- -1.0
for(i in seq_len(nrow(df_count_table))) {
  mat <- df_count_table[i,c('call_found','claim_found','call_notfound','claim_notfound')] %>% 
    as.numeric() %>%
    matrix(nrow=2, byrow=T)
  # res <- fisher.test(mat, alternative='g')
  res <- fisher.test(mat) # two.sided
  df_count_table[i,]$pvalue <- res$p.value
}

# One-sided: n=228 (p<0.05), n=40 (p>=0.05)
# Two-sided: n=235 (p<0.05), n=33 (p>=0.05)
dplyr::count(df_count_table, pvalue < 0.05)
any(df_count_table$pvalue == -1.0)

subset(df_count_table,pvalue > 0.05)
subset(df_count_table,pvalue < 1e-5)



# Pure contamination (homogenous) - main boxplot
df_var_misid <- subset(df_var, order %in% samples_misid)
length(samples_misid) # 3,046
length(samples_cross) #  248

df_var_misid_median <- df_var_misid %>%
  group_by(sample, group) %>%
  summarize(median=median(vaf_obs))

median(subset(df_var_misid_median, group=='call - 1')$median) # 0.53395
median(subset(df_var_misid_median, group=='call - 2')$median) # 1.0


# Figure. CCLE2CCLE VAF by boxplot, misid only (n=3,046). claim nSNP < 3. median = 0.53395 (hetero), 1.0 (homo)
set.seed(1)
center_stroke <- df_var_misid_median %>%
  group_by(group) %>%
  summarize(median=quantile(median,0.5))
mylabels2 <- c('call - 1'= 'Hetero (Call)', 'call - 2'= 'Homo (Call)',
              'claim - 1'= 'Hetero (Claim)', 'claim - 2'= 'Homo (Claim)')
mylabels2 <- factor(mylabels2, levels=c('Hetero (Call)','Homo (Call)','Hetero (Claim)','Homo (Claim)'))

ggplot(df_var_misid_median, aes(x=mylabels2[group], y=median, fill=group)) +
  geom_violin(color=NA, show.legend=F, alpha=1, lwd=1) +
  geom_boxplot(data=center_stroke, stat='identity', # Median mark white
               aes(ymin=median, ymax=median, # No wiskers
                   lower=median, middle=median, upper=median), 
               color='white', fill=NA, lwd=0.5, show.legend = F) +
  geom_jitter(color='black', size=0.1, alpha=0.1, show.legend=F) +
  scale_fill_manual(values=mypalette, labels=mylabels, name='genotype') +
  coord_cartesian(ylim=c(0.0,1)) +
  xlab('') +
  ylab('median VAF') +
  theme_pubr() +
  theme(legend.position='right', 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename=paste0(path_root, '/pdfs/20241111_FigureX_ccle2ccle_vaf_by_boxplot_misid_both.pdf'), unit="in",width=3, height=4)




# Pure contamination (homogenous) - supplementary
df_var_misid <- subset(df_var, order %in% samples_misid & cell == 'call')
df_var_misid_mixture <- unique(df_var_misid[,c('order','cell','frac')])
df_var_misid_stat <- df_var_misid %>%
  group_by(order, group, cell) %>%
  summarize(ymin = quantile(vaf_obs, 0.00),
            lower = quantile(vaf_obs, 0.25),
            median = median(vaf_obs),
            upper = quantile(vaf_obs, 0.75),
            ymax = quantile(vaf_obs, 1.00),
            total=n())
df_misid_count <- df_var_misid %>%
  dplyr::count(order, cell)

draw_misid <- function(sample_list) { # call or claim
  # subset from global variables
  tmp_df_var <- subset(df_var_misid, order %in% sample_list)
  tmp_df_var_stat <- subset(df_var_misid_stat, order %in% sample_list)
  tmp_df_var_mixture <- subset(df_var_misid_mixture, order %in% sample_list)
  tmp_df_count <- subset(df_misid_count, order %in% sample_list)
  
  plot_base <- ggplot(tmp_df_var) +
    scale_fill_manual(values=mypalette, labels=mylabels, name='genotype') +
    scale_color_manual(values=mypalette, labels=mylabels, name='genotype') +
    theme_pubr() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),)
  
  plot_VAF <- plot_base +
    geom_hline(yintercept=c(0, 0.5, 1), linetype=2, alpha=0.25) +
    geom_boxplot(data=tmp_df_var,
                 outlier.shape = NA, position = position_dodge2(preserve = "single"), lwd=0.5,
                 aes(x=order, y=vaf_obs, group = interaction(order, group),
                     fill=group, color=group), alpha=1.0) +
    geom_boxplot(data=tmp_df_var_stat, stat='identity',
                 outlier.shape = NA, position = position_dodge2(preserve = "single"), # This is for median
                 aes(ymin=median, ymax=median, # No wiskers
                     x=order, lower=median, middle=median, upper=median, group = interaction(order, group),
                 ), color='white', fill=NA, lwd=0.5, show.legend = F) +
    geom_line(data=tmp_df_var_mixture, aes(x = order, y = frac, group=1), color='red', linetype=1) +
    ylab('Varient allele frequency') + xlab('') +
    theme(axis.title.x=element_blank())
  
  plot_nSNP <- plot_base + 
    geom_hline(yintercept=c(10,100,1000), linetype=2, alpha=0.25) +
    geom_bar(data=tmp_df_count, stat='identity', width=1,
             aes(x=order, y=n)) +
    scale_y_log10(limits=c(1,5000)) +
    ylab('nSNP') + xlab('Samples')
  
  return (ggarrange(plot_VAF, plot_nSNP, ncol=1, nrow=2, align='v', heights=c(4,2.5), common.legend = T, legend='right'))
}


# Figure. CCLE2CCLE VAF by boxplot, mis-identification only (n=3,046). claim nSNP < 3. - order (-frac). call shown. 500 per each row. 1-500, 501-1000, ..., 2500-3046
# each 500 records
samples_misid_sort <- samples_misid[order(samples_misid)] # sort by level
bin <- c(0, 500, 1000, 1500, 2000, 2500, 3046)
plots <- list()
for (i in 1:(length(bin)-1)) {
  start = bin[i] + 1
  end = bin[i+1]
  samples_misid_subset <- samples_misid_sort[start:end]
  plots[[i]] <- draw_misid(samples_misid_subset)
}
ggarrange(plotlist=plots, ncol=1, nrow=length(bin)-1, align='v')

ggsave(filename=paste0(path_root, '/pdfs/20240910_FigureX_ccle2ccle_vaf_by_boxplot_misid_supplementary_500.pdf'), unit="in",width=12.0, height=18.0)

