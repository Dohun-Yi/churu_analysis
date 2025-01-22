# contamination_profile
# Responsible for...
#   contamination profile for US, China, others

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

path_root <- "../../"
# using representative name for visualize purpose
path_country <- paste0(path_root, '/analysis/3_country_institute/institute_selected.parsed.corrected.add_time.txt.gz') # latest
path_contam <- paste0(path_root, '/analysis/4_check_misidentified/1697001520.compare_result.final.representative.txt') # latest
path_institute_agg <- paste0(path_root, '/analysis/4_check_misidentified/misidentified_institute.txt') 
path_journal <- paste0(path_root, '/analysis/5_journals/sam2journal.txt') # latest

df_contam <- fread(path_contam, header=F) %>% as.data.frame()
df_country <- fread(path_country, header=F) %>% as.data.frame()
df_institute_agg <- fread(path_institute_agg, header=F) %>% as.data.frame()
df_journal <- fread(path_journal, header=F) %>% as.data.frame()
df_merge <- merge(df_contam, df_country, by='V1', all.x=T)
df_merge <- merge(df_merge, df_journal, by='V1', all.x=T)
df_merge <- merge(df_merge, df_institute_agg, by.x='V3.y', by.y='V1', all.x=T)

names(df_merge) <- c('RorID', 'Sample','Call','Claim','-','Match_raw',
                     'Owner', '_1', '_2', 'Relation', 'Submission',
                     "Project", "PMID", "Journals", "IFs", "Journal", "IF",
                     "Institute", "Country", "_3", "_4", "_5", "Frac_institute",
                     "_6", "_7", "_8", "Frac_aggregation", 'Parents')

df_merge$Match <- str_remove_all(df_merge$Match_raw,'-.*')
df_merge[df_merge$Match == 'contam',]$Match <- 'mismatch' # for convenience
df_merge_filt <- subset(df_merge, Match != "skip" & !is.na(Country))

nrow(df_merge_filt) # 79,638 (records without contries are not required here)

df_merge_filt$Claim <- str_replace_all(df_merge_filt$Claim, '_primary','primary') # just renaming
# Note that, this specific analysis include nonCCLE contaminations, since there's no percentage here

# ============================================================== #

df_country_size <- dplyr::count(df_merge_filt,Country)
df_country_size <- df_country_size[order(df_country_size$n,decreasing=T),]
top4 <- df_country_size[1:4,]$Country # Countries with top 4 data submission = US, China, UK, Germany

df_subset_mismatch <- subset(df_merge_filt, Match=='mismatch')[,c('Call','Claim','Country', 'Institute')]
df_profile <- dplyr::count(df_subset_mismatch, Call, Claim, Country, Institute)
df_profile <- df_profile[order(df_profile$n, decreasing=T),]
df_profile$pattern <- paste(sep=' > ', df_profile$Claim, df_profile$Call)


draw_country <- function(df, country, top_n=20) {
  if (country=='All others') {
    df_tmp <- df
  } else {
    df_tmp <- subset(df, Country == country)
  }
  
  # Find top pattern
  df_sum <- df_tmp %>% 
    group_by(pattern) %>% 
    summarize(sum=sum(n))
  top_pattern <- df_sum[order(df_sum$sum,decreasing=T),]$pattern[1:top_n]
  
  # subset, sort, and add color group
  df_tmp <- subset(df_tmp, pattern %in% top_pattern)
  df_tmp$pattern <- factor(df_tmp$pattern, levels=rev(top_pattern))
  df_tmp$group <- paste0(df_tmp$pattern,'-',df_tmp$Institute)
  df_tmp$group <- factor(df_tmp$group, levels=rev(df_tmp$group))
  
  # make color palette
  palette_size <- 30 # maximum is 30 (China, Hep-G2 > HeLa)
  fixed_palette <- gradient_n_pal(c("#BBBBBB", "#F3CEBB", "#B31312"), c(0,5,30))
  
  patterns <- df_tmp$pattern
  color_scale_list <- list()
  for (i in seq_along(patterns)){
     df_sub <- subset(df_tmp, pattern==patterns[i])
     # df_sub$group <- factor(df_sub$group, levels=df_sub$group)
     color_scale <- fixed_palette(1:nrow(df_sub))
     names(color_scale) <- df_sub$group
     color_scale_list <- append(color_scale_list, color_scale)
  }
  ggplot(df_tmp, aes(y=pattern, x=n, fill=group)) +
    geom_bar(stat='identity', position='stack', color='black', size=1, width=1, show.legend=F) +
    ylab('') + xlab('Number of samples') + ggtitle(country) +
    scale_fill_manual(values=color_scale_list) +
    scale_x_continuous(breaks = breaks_pretty()) +
    theme_pubr()
}

# Figure. Contamination profile for countries with top two data submission. segmented by different institute, and color is number of institute.
p1 <- draw_country(df_profile, 'United States', 20)
p2 <- draw_country(df_profile, 'China', 20)
p3 <- draw_country(subset(df_profile, ! Country %in% c('United States','China')), 'All others', 20)
ggarrange(p1, p2, p3, align='h', ncol=3, nrow=1)
ggsave(filename=paste0(path_root, '/pdfs/20241103_FigureX_contamination_profile_UCO.pdf'), unit="in", width=12.0, height=4.0)
