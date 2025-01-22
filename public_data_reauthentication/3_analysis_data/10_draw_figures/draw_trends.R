# draw_trends
# Responsible for...
#   contamination in journals
#   yearly contamination

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
path_country <- paste0(path_root, '/analysis/3_country_institute/institute_selected.parsed.corrected.add_time.txt.gz') # latest
path_contam <- paste0(path_root, '/analysis/4_check_misidentified/1697001520.compare_result.final.txt') # latest
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

df_merge_filt <- subset(df_merge, Match != "skip")
nrow(df_merge_filt) # 79,696 using all data. regardless of country is omitted or not

df_merge_filt <- subset(df_merge_filt, ! Match_raw %in% c("contam-unknown-nonccle", "contam-known-nonccle"))
nrow(df_merge_filt) # 78,575 excluding all nonCCLE cell lines



# ============================================================== #
# Analysis: Recent trends in China 20241106
df_merge_filt_china <- subset(df_merge_filt, Country=='China')
df_merge_filt_china$year <- format(df_merge_filt_china$Submission, '%Y')
df_merge_filt_china <- subset(df_merge_filt_china, !is.na(year))
df_merge_filt_china$year <- format(df_merge_filt_china$Submission, '%Y')
df_merge_filt_china[df_merge_filt_china$year %in% c('2010','2011','2012'),]$year <- '-2012'
df_merge_filt_china[df_merge_filt_china$year %in% c('2022','2023'),]$year <- '2022-'

df_frac <- df_merge_filt_china %>% 
  group_by(year) %>%
  summarize(n_mat=sum(Match=='match'),
            n_mis=sum(Match=="mismatch"))
df_frac$total <- df_frac$n_mat + df_frac$n_mis
df_frac$frac <- df_frac$n_mis / df_frac$total

# show table
show(df_frac)



# ============================================================== #
# Analysis: Recent trend, all countries

# Need data submission year information
df_merge_filt$year <- format(df_merge_filt$Submission, '%Y')
df_merge_filt$yearbin <- ifelse(as.numeric(df_merge_filt$year) < 2017, 'before 2017', 'since 2017')

# Some data is omitted on "meta" file. need correction later (got removed in MERGE step)
# nrow(df_merge_filt) # now 79696 because "unsure" is removed
nrow(df_merge_filt) # now 78575 because nonCCLE is further removed

# Filtration
df_count_tmp <- dplyr::count(df_merge_filt,Country)
countries_of_interest <- subset(df_count_tmp,n >= 100)$Country
df_merge_filt_interest <- subset(df_merge_filt, Country != 'DISCARD' & Country %in% countries_of_interest)

df_frac <- df_merge_filt_interest %>% 
  group_by(Country, yearbin) %>%
  summarize(n_mat=sum(Match=='match'),
            n_mis=sum(Match=="mismatch"))
df_frac$total <- df_frac$n_mat + df_frac$n_mis
df_frac$frac <- df_frac$n_mis / df_frac$total
df_frac <- df_frac[order(df_frac$frac, decreasing = T),] # sort

df_frac[df_frac$Country=='United States of America',]$Country <- 'United States' # Just for visualization
country_order_tmp <- subset(df_frac, yearbin=='since 2017')
country_order_tmp <- country_order_tmp[order(country_order_tmp$frac, decreasing=T),]$Country
df_frac$Country <- factor(df_frac$Country, levels=country_order_tmp)
df_frac <- df_frac %>% as.data.frame()
df_frac$legend <- '<2017'
df_frac[df_frac$yearbin=='since 2017',]$legend <- '>=2017'


# Figure. contamination trends with data size - Final (hopefully)
ggplot(df_frac, aes(y=frac*100, x=Country, color=factor(legend))) +
  geom_line(aes(group=Country), color='black') +
  geom_point(data=subset(df_frac, legend=='>=2017'), shape=19, aes(size=log10(total))) +
  geom_point(data=subset(df_frac, legend=='<2017'), shape=19, aes(size=log10(total))) +
  theme_pubr() +
  xlab('') + ylab('Percentage of misidentified data (%)') +
  ggtitle('')+
  scale_color_manual(values=c('<2017'='#26252C','>=2017'='#E54861'), name='Year') +
  scale_size_continuous(breaks=c(1,2,3,4), labels=c(10,100,1000,10000),name='Data size', 
                        limits = c(1,6), range=c(1,7)) +
  theme(legend.position = 'right', axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename=paste0(path_root, '/pdfs/20241024_FigureX_contamination_trends_size.pdf'), unit="in",width=9.0, height=4.0)




# ============================================================== #
# Analysis: Recent trend, overall

df_merge_filt_yealy <- subset(df_merge_filt, !is.na(year))
df_merge_filt_yealy$year <- format(df_merge_filt_yealy$Submission, '%Y')
df_merge_filt_yealy[df_merge_filt_yealy$year %in% c('2010','2011','2012'),]$year <- '-2012'
df_merge_filt_yealy[df_merge_filt_yealy$year %in% c('2022','2023'),]$year <- '2022-'
count(subset(df_merge_filt_yealy, year %in% c('-2012','2013','2014','2015','2016')), Match)
count(subset(df_merge_filt_yealy, year %in% c('2017','2018','2019','2020','2021','2022-')), Match)

df_merge_filt_yealy$Label <- "Correct"
df_merge_filt_yealy[df_merge_filt_yealy$Match=="mismatch",]$Label <- "Misidentified"

# Figure. Yearly data submission and contamination (year binned) - Final (hopefully)
p1 <- ggplot(df_merge_filt_yealy, aes(fill=Label, x=year)) +
  geom_bar(stat='count', position='stack') +
  scale_fill_manual(values=c('Correct'='#26252C','Misidentified'='#E54861'), name='') +
  ylab('Number of samples') + xlab('') +
  theme_pubr() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- ggplot(df_merge_filt_yealy, aes(fill=Label, x=year)) +
  geom_bar(stat='count', position='fill') +
  scale_fill_manual(values=c('Correct'='#26252C','Misidentified'='#E54861'), name='') +
  scale_y_continuous(labels=scales::percent, breaks=c(0:5*0.02)) +
  coord_cartesian(ylim=c(0,0.10)) +
  xlab('Year') + ylab("Percentage of misidentified data") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggarrange(p1, p2, align='v', nrow=2, ncol=1, common.legend = T, legend = 'right', heights = c(1.60,2))
ggsave(filename=paste0(path_root, '/pdfs/20241024_FigureX_yearly.pdf'), unit="in",width=6.0, height=6.0)



# ============================================================== #
# Analysis: Journal
MINIMAL_DATA <- 100

df_merge_filt$year <- format(df_merge_filt$Submission, '%Y')
df_merge_filt_journal <- subset(df_merge_filt, ! is.na(IF))
df_merge_filt_journal$journal <- df_merge_filt_journal$Journal
df_merge_filt_journal$ifbin <- as.integer(df_merge_filt_journal$IF/10)*10 # impact factor
df_merge_filt_journal[df_merge_filt_journal$ifbin > 50, ]$ifbin <- 50 # impact factor
length(unique(df_merge_filt_journal$Journal)) 
nrow(df_merge_filt_journal)

df_frac <- df_merge_filt_journal %>% 
  group_by(journal, IF) %>%
  summarize(n_mat=sum(Match=='match'),
            n_mis=sum(Match=="mismatch"))
df_frac$total <- df_frac$n_mat + df_frac$n_mis
df_frac$frac <- df_frac$n_mis / df_frac$total
df_frac <- df_frac[order(df_frac$frac, decreasing = T),] # sort
df_frac[df_frac$journal=="Scientific data",]$IF <- 9.8 # Correction for false high-impact (?)
df_frac[df_frac$journal=="Science (New York, N.Y.)",]$journal <- "Science"
df_frac_all <- df_frac
df_frac_filt <- df_frac[df_frac$total >= MINIMAL_DATA,] # 



# Figure. Misidentified data by journals, more than 100 data
ggplot(df_frac, aes(x=IF, y=frac*100, size=log10(total), color=log10(total))) +
  geom_point(data=subset(df_frac, journal != "" & total>=1000)) +
  geom_point(data=subset(df_frac, journal != "" & total< 1000 & total > 100)) +
  scale_x_continuous(breaks=c(0:100*10), limits=(c(0,65))) +
  scale_y_continuous(breaks=c(0:10*5), limits=(c(0,30))) +
  scale_size_continuous(name = 'Data size',
                        breaks = 2:4, labels = 10**c(2:4),
                        limits = c(2,4), range=c(1.5,5)) +
  scale_color_gradientn(name = 'Data size',
                        colors = c('#000000','#E6AC00','#A52A2A'),
                        breaks = 2:4, labels = 10**c(2:4),
                        limits = c(2,4)) +
  xlab('Impact factor (3 years)') +
  ylab('Percentage of misidentified data (%)') +
  theme_classic2() +
  theme(panel.grid.major = element_line(colour = "grey", linetype=2)) +
  guides(color = guide_legend(ncol = 1), size = guide_legend(ncol = 1)) +
  theme(legend.position = c(0.7,0.7),
        legend.background = element_rect(fill = "white", color = "black"))

ggsave(filename=paste0(path_root, '/pdfs/20241024_FigureX_contamination_journal.pdf'), unit="in",width=5.0, height=3.5)



# Top three journals
df_meta_filt_CNS <- subset(df_merge_filt_journal, journal %in% c('Nature', 'Cell' ,'Science (New York, N.Y.)'))
df_meta_filt_CNS[1,]
df_meta_filt_CNS$journal <- df_meta_filt_CNS$Journal
df_frac <- df_meta_filt_CNS %>% 
  group_by(journal, PMID, year) %>%
  summarize(n_mat=sum(Match=='match'),
            n_mis=sum(Match=="mismatch"))
df_frac$total <- df_frac$n_mat + df_frac$n_mis
df_frac$frac <- df_frac$n_mis / df_frac$total
df_frac <- df_frac[order(df_frac$frac, decreasing = T),] # sort
for (i in 2012:2023) {
  for (j in unique(df_frac$journal)) {
    if (nrow(subset(df_frac, journal==j & year==i))==0){
      df_frac = rbind(df_frac, data.frame(
        journal=j,
        year=as.character(i),
        n_mat=0,
        n_mis=0,
        total=0,
        frac=0
      ))
    }
  }
}
df_frac_m <- melt(df_frac, measure.vars=c('n_mis','n_mat'))
df_frac_m[grepl("Science",df_frac_m$journal),]$journal <- "Journal A"
df_frac_m[grepl("Nature",df_frac_m$journal),]$journal <- "Journal B"
df_frac_m[grepl("Cell",df_frac_m$journal),]$journal <- "Journal C"
df_frac_m$var <- "Correct"
df_frac_m[grepl("n_mis",df_frac_m$variable),]$var <- "Misidentified"
df_frac_m$var <- factor(df_frac_m$var, levels=c('Correct','Misidentified'))


# Figure. Misidentified data by top three journals (CNS) - Final (hopefully)
ggplot(subset(df_frac_m, year!=2011 & var=="Misidentified"), aes(x=year, y=value, fill=journal)) +
  geom_bar(stat='identity', position='stack', color='white') +
  scale_fill_manual(values=c("Journal A"='#E59500',
                             "Journal B"='#840032',
                             "Journal C"='#002642'
  ), name='') +
  xlab('Year of publication') +
  ylab('Number of misidentified data') +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1),
        legend.position = 'right')

ggsave(filename=paste0(path_root, '/pdfs/20241024_FigureX_contamination_journal_top3.pdf'), unit="in",width=6.0, height=3.0)



# Checking authentication rates
# Journals and its authentication rates are collected from (https://doi.org/10.1016/j.isci.2020.101698)

path_auth_rate <- paste0(path_root, 'analysis/5_journals/iscience_suptable.csv')
df_auth_rate <- fread(path_auth_rate, header=T) %>% as.data.frame()
names(df_auth_rate) <- c('paper','year','journal','v1','v2','auth')
df_auth_rate$authenticated <- round(df_auth_rate[,1] * as.numeric(str_remove_all(df_auth_rate[,6],'%'))/100)

df_auth_total <- df_auth_rate %>%
  group_by(journal) %>%
  summarize(total=sum(paper),
            auth=sum(authenticated))
df_auth_yearly <- df_auth_rate %>%
  group_by(year) %>%
  summarize(total=sum(paper),
            auth=sum(authenticated))
df_auth_total$authrate <- df_auth_total$auth / df_auth_total$total
df_auth_yearly$authrate <- df_auth_yearly$auth / df_auth_yearly$total

df_yearly_journal <- df_merge_filt_journal %>% 
  subset(journal %in% df_auth_rate$journal) %>%
  group_by(year) %>%
  summarize(n_mat=sum(Match=='match'),
            n_mis=sum(Match=="mismatch"),
            frac=n_mis/(n_mat+n_mis))

length(unique(df_frac_filt$journal)) # 79 journals with 100 data submission
length(unique(df_auth_total$journal)) # 79 journals with 100 data submission
length(df_frac_filt[! df_frac_filt$journal %in% unique(df_auth_rate$Title),]$journal) # 18 unavailable

df_merge <- merge(df_frac_filt,df_auth_total, by='journal')
nrow(df_merge) # 61 or 206

# Figure. authentication rate and contamination rate.
df_merge_journal_year <- merge(df_auth_yearly, df_yearly_journal, by='year')
sum(df_merge_journal_year$n_mat + df_merge_journal_year$n_mis) # 20,823

ggplot(df_merge_journal_year, aes(x=year)) +
  geom_point(aes(y=authrate*100, color='paper')) +
  geom_line(aes(y=authrate*100, color='paper')) +
  geom_point(aes(y=frac*100, color='data'), ) +
  geom_line(aes(y=frac*100, color='data')) +
  scale_color_manual(values=c("paper"='#26252C',"data"='#E54861'),
                     labels=c("paper"='Authenticated papers (%)',"data"='Contaminated data (%)'),
                     name='') +
  xlab('Year') + ylab('Percentage') +
  ylim(c(0,15)) +
  scale_x_continuous(breaks=c(2010:2019)) +
  theme_pubr() +
  theme(legend.position = c(0.3,0.8), legend.background = element_blank(),
        axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))

ggsave(filename=paste0(path_root, '/pdfs/20241101_FigureX_auth_and_contam.pdf'), unit="in",width=5.0, height=3.5)



