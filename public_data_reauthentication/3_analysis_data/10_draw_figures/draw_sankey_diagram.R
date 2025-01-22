# sankey_diagram
# Responsible for...
#   Sankey diagram (contamination patterns)

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
library(ggsankey)
library(pals)

path_root <- "../../"
path_country <- paste0(path_root, '/analysis/3_country_institute/institute_selected.parsed.corrected.add_time.txt.gz') # latest
path_contam <- paste0(path_root, '/analysis/4_check_misidentified/1697001520.compare_result.final.representative.txt') # latest

df_contam <- fread(path_contam, header=F) %>% as.data.frame()
df_country <- fread(path_country, header=F) %>% as.data.frame()
df_merge <- merge(df_contam, df_country, by='V1', all.x=T)

names(df_merge) <- c('Sample','Call','Claim','-','Match_raw',
                     'Owner', '_1', '_2', 'Country', 'Relation', 'Submission')

df_merge$Match <- str_remove_all(df_merge$Match_raw,'-.*')
df_merge[df_merge$Match == 'contam',]$Match <- 'mismatch' # for convenience
df_merge_filt <- subset(df_merge, Match != "skip" & !is.na(Country))

nrow(df_merge_filt) # 79,641 (records without contries are not required here)

df_merge_filt$Claim <- str_replace_all(df_merge_filt$Claim, '_primary','primary') # just renaming


# Draw
draw_by_country <- function(df, countries, title, n_sample=20, reverse = F) {
  if (!reverse) {
    df <- subset(df, Country %in% countries)[,c('Call','Claim')]
  } else {
    df <- subset(df, ! Country %in% countries)[,c('Call','Claim')]
  }
  
  names(df) <- c('target','source')
  
  target_high <- (dplyr::count(df, target) %>% slice_max(n, n=n_sample))$target
  source_high <- (dplyr::count(df, source) %>% slice_max(n, n=n_sample))$source
  
  if (any(! df$target %in% target_high)) {
    df[! df$target %in% target_high,]$target <- 'others'
  }
  if (any(! df$source %in% source_high)) {
    df[! df$source %in% source_high,]$source <- 'others'
  }
  
  df_long <- df %>%
    make_long(source, target)
  
  name_order <- names(sort(table(c(df$target, df$source))))
  name_order <- c("others",name_order[1:length(name_order)-1])
  df_long$node <- factor(df_long$node, levels=name_order)
  df_long$next_node <- factor(df_long$next_node, levels=name_order)
  
  ncolors <- length(name_order)
  ggplot(df_long, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
    geom_sankey(space=50, flow.alpha=0.5) +
    geom_sankey_label(space=50, size = 3, fill='white') +
    scale_fill_manual(values = ocean.thermal(n=ncolors+3)[1:ncolors]) +
    ggtitle(title) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
}

df_origin <- subset(df_merge_filt,Match=='mismatch')[,c('Call','Claim', 'Country')]

n_sample <- 10
p1 <- draw_by_country(df_origin, 'United States', 'United States', n_sample=n_sample)
p2 <- draw_by_country(df_origin, 'China', 'China', n_sample=n_sample)
p3 <- draw_by_country(df_origin, c('United States', 'China'), 'Others', n_sample=n_sample, reverse=T)

# Figure. Sankey for US, China, and others
ggarrange(p1, p2, p3, ncol=3, nrow=1, align='hv')
ggsave(filename=paste0(path_root, '/pdfs/20241103_FigureX_sankey_10.pdf'), unit="in",width=20.0, height=8.0)
