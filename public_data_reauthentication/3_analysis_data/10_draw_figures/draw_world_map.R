# world_map
# Responsible for...
#   Contamination map
#   Data acquisition map

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
library(factoextra)
library(ggrepel)
library(tools)
library(gridExtra)

library(sf)
library(rmapshaper)
library(rworldmap)
library(rnaturalearth)

library(grid)
library(rsvg)
library(magick)

# Data load
path_root <- "../../"
path_country <- paste0(path_root, '/analysis/3_country_institute/institute_selected.parsed.corrected.add_time.txt.gz') # latest
path_contam <- paste0(path_root, '/analysis/4_check_misidentified/1697001520.compare_result.final.txt') # latest
path_institute_agg <- paste0(path_root, '/analysis/4_check_misidentified/misidentified_institute.txt') 

df_contam <- fread(path_contam, header=F) %>% as.data.frame()
df_country <- fread(path_country, header=F) %>% as.data.frame()
df_institute_agg <- fread(path_institute_agg, header=F) %>% as.data.frame()
df_merge <- merge(df_contam, df_country, by='V1', all.x=T)
df_merge <- merge(df_merge, df_institute_agg, by.x='V3.y', by.y='V1', all.x=T)
names(df_merge) <- c('RorID', 'Sample','Call','Claim','-','Match_raw',
                     'Owner', '_1', '_2', 'Relation', 'Submission',
                     "Institute", "Country", "_3", "_4", "_5", "Frac_institute",
                     "_6", "_7", "_8", "Frac_aggregation", 'Parents')
df_merge$Match <- str_remove_all(df_merge$Match_raw,'-.*')
df_merge[df_merge$Match == 'contam',]$Match <- 'mismatch' # for convenience

df_merge_filt <- subset(df_merge, Match != "skip" & !is.na(Country))
nrow(df_merge_filt) # 79,638 (records without contries are not required here) - Was 79,641, but reduced to 79,638 at some time point

df_merge_filt <- subset(df_merge_filt, ! Match_raw %in% c("contam-unknown-nonccle", "contam-known-nonccle"))
nrow(df_merge_filt) # 78,520 excluding all nonCCLE cell lines (updated 202410)

nrow(subset(df_merge, Match != "skip" & is.na(Country))) # 58 records are removed because Country is NA (due to error or something)

# ============================================================== #
# Refine country name - Natural earth naming convention
df_merge_filt[df_merge_filt$Country=='United States',]$Country <- 'United States of America'


# Let's check some "DISCARD"
df_discard <- subset(df_merge_filt, Match=='mismatch' & Country == "DISCARD")
df_discard_count <- dplyr::count(df_discard,Owner)
df_discard_count <- df_discard_count[order(df_discard_count$n, decreasing = T),] # sort
show(df_discard_count)


# Data processing
MINIMAL_DATA <- 100 # minimal data number=100


# Get contamination fraction
df_frac <- df_merge_filt %>% 
  group_by(Country) %>%
  summarize(n_mat=sum(Match=='match'),
            n_mis=sum(Match=="mismatch"))
df_frac$total <- df_frac$n_mat + df_frac$n_mis
df_frac$frac <- df_frac$n_mis / df_frac$total
df_frac <- df_frac[order(df_frac$frac, decreasing = T),] # sort
df_frac_filt <- df_frac[df_frac$total >= MINIMAL_DATA,] # 


# World map data
world <- ne_countries(scale = 10, returnclass = "sf")
dispute <- ne_download(scale = 10, type = "admin_0_disputed_areas", category = "cultural")

world <- world %>% filter(iso_a3 != "ATA") # remove antarctica
world <- world %>%
  left_join(df_frac_filt, by = c("name" = "Country"))

world_simplified <- ms_simplify(world, keep = 0.01, keep_shapes = TRUE)
dispute_simplified <- ms_simplify(dispute, keep = 0.01, keep_shapes = TRUE)

robinson_proj <- "+proj=robin +lon_0=0 +datum=WGS84"
world_robinson <- st_transform(world_simplified, crs = robinson_proj)
dispute_robinson <- st_transform(dispute_simplified, crs = robinson_proj)

subset(df_frac_filt, Country %in% c('Singapore', 'Israel'))$frac
small_countries <- data.frame(
  name = c("Singapore", "Israel"),
  lon = c(103.8198, 34.8516),
  lat = c(1.3521, 31.0461),
  frac = c(0.09178082, 0.07196970) # Need to be changed when data updated!! (last update 20241008)
)

# Transform small country coordinates to Robinson projection
small_countries_sf <- st_as_sf(small_countries, coords = c("lon", "lat"), crs = 4326)
small_countries_robinson <- st_transform(small_countries_sf, crs = robinson_proj)

# Is there any missing country name?
df_frac_filt$Country[! df_frac_filt$Country %in% world$name]

palette <- c('#666666', '#FFD23F', '#9B0A0A')
value_scale <- c(0.0, 0.25, 0.5, 0.75, 1.0)


# 20240907. Figure. Contamination map - Europe inset plot added (10 x 4 inch) - Final (hopefully)
plot_europe <- ggplot(data = world_robinson) +
  geom_sf(aes(fill = frac*100), color='black', show.legend = F) + # Contamination
  geom_sf(data=subset(world_robinson, is.na(frac)),fill='#EEEEEE', color='black') + # Empty points
  geom_sf(data=dispute_robinson, fill='#BBBBBB', color=NA, linetype=2) + # Border dispute
  scale_fill_gradientn(colors = palette, name='Contaminated (%)', 
                       values=value_scale) + 
  xlim(c(-1000000, 3000000)) +
  ylim(c(3600000, 7400000)) +
  theme_void() +
  theme(panel.background = element_rect(fill='white', color='black'))

ggplot(data = world_robinson) +
  geom_sf(aes(fill = frac*100), color='black') + # Contamination
  geom_sf(data=subset(world_robinson, is.na(frac)),fill='#EEEEEE', color='black') + # Empty points
  geom_sf(data=dispute_robinson, fill='#BBBBBB', color=NA, linetype=2) + # Border dispute
  geom_sf(data=small_countries_robinson, aes(fill=frac*100), color='black', shape=21, size=3) + # Small countries
  scale_fill_gradientn(colors = palette, name='Contaminated (%)', 
                       values=value_scale) +
  annotation_custom(
    grob = ggplotGrob(plot_europe),
    xmin = -4000000,
    xmax = +5000000,
    ymin = -7000000, 
    ymax = +3000000
  ) +
  annotate("rect",fill=NA, color='black', linewidth=0.5,
           xmin=-1000000, xmax=3000000,
           ymin=3600000, ymax=7500000)+
  annotate("rect",fill=NA, color='black', linewidth=0.5,
           xmin=-4000000, xmax=5000000,
           ymin=-6270000, ymax=2260000)+
  annotate("segment", colour = "black", linetype=2,
           x = -4000000, xend = -1000000, 
           y = 2260000, yend = 3600000) + # Left
  annotate("segment", colour = "black", linetype=2,
           x = 3000000, xend = 5000000, 
           y = 3600000, yend = 2260000) + # Right
  theme_void() +
  guides(fill = guide_colorbar(ticks.colour = "white")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.95,0.7))
  
ggsave(filename=paste0(path_root, '/pdfs/20241111_FigureX_contamination_map.pdf'), unit="in",width=10.0, height=4.0)



# ======== Data Collection Map ======== 

# Data load
path_root <- "../../"
path_ror <- paste0(path_root, '/analysis/3_country_institute/table') # latest
path_country <- paste0(path_root, '/analysis/3_country_institute/institute_selected.parsed.corrected.add_time.txt.gz') # latest
path_contam <- paste0(path_root, '/analysis/4_check_misidentified/1697001520.compare_result.final.txt') # latest
path_institute_agg <- paste0(path_root, '/analysis/4_check_misidentified/misidentified_institute.txt') 


df_contam <- fread(path_contam, header=F) %>% as.data.frame()
df_country <- fread(path_country, header=F) %>% as.data.frame()
df_institute_agg <- fread(path_institute_agg, header=F) %>% as.data.frame()
df_ror <- fread(path_ror, header=T) %>% as.data.frame()
df_coordinate <- df_ror[,c('id', 'addresses[0].lat', 'addresses[0].lng')]
df_merge <- merge(df_contam, df_country, by='V1')
df_merge <- merge(df_merge, df_institute_agg, by.x='V3.y', by.y='V1')
df_merge_coord <- merge(df_merge, df_coordinate, by.x='V3.y', by.y='id', all.x=T)
names(df_merge_coord) <- c('ROR', 'Sample','Call','Claim','-','Match_raw',
                           'Owner', '_1', '_2', 'Relation', 'Submission',
                           "Institute", "Country", "_3", "_4", "_5", "Frac_institute",
                           "_6", "_7", "_8", "Frac_aggregation","_9",
                           'lat', 'lon')

df_merge_coord$Match <- str_remove_all(df_merge_coord$Match_raw,'-.*')
df_merge_coord[df_merge_coord$Match == 'contam',]$Match <- 'mismatch' # for convenience
df_merge_coord <- subset(df_merge_coord, Match != "skip")
df_merge_filt <- subset(df_merge_coord, Match %in% c('match', 'mismatch'))
nrow(df_merge_filt) # 79,638. This parts includes nonCCLE

# ============================================================== #
# Refine country name - Natural earth naming convention
df_merge_filt[df_merge_filt$Country=='United States',]$Country <- 'United States of America'


# Get total data size
df_total <- dplyr::count(df_merge_filt, Country)


# Get total number of samples per coordinate - some institutes share coordination, so just coord level aggregation
df_coord_count <- dplyr::count(df_merge_filt, lat, lon) %>%
  subset(!is.na(lat))


# World map data
world <- ne_countries(scale = 10, returnclass = "sf")
dispute <- ne_download(scale = 10, type = "admin_0_disputed_areas", category = "cultural")

world <- world %>% filter(iso_a3 != "ATA") # remove antarctica
world <- world %>%
  left_join(df_total, by = c("name" = "Country"))

world_simplified <- ms_simplify(world, keep = 0.01, keep_shapes = TRUE)
dispute_simplified <- ms_simplify(dispute, keep = 0.01, keep_shapes = TRUE)

robinson_proj <- "+proj=robin +lon_0=0 +datum=WGS84"
world_robinson <- st_transform(world_simplified, crs = robinson_proj)
dispute_robinson <- st_transform(dispute_simplified, crs = robinson_proj)

df_coord_count_sf <- st_as_sf(df_coord_count, coords = c("lon", "lat"), crs = 4326)
df_coord_count_robinson <- st_transform(df_coord_count_sf, crs = robinson_proj)

coord_shift <- function(df, shift) {
  df_shift <- df
  st_geometry(df_shift) <- st_geometry(df_shift) + c(0,shift)
  st_crs(df_shift) <- robinson_proj
  df_shift
}

df_shift1 <- coord_shift(df_coord_count_robinson, 52000)
df_shift2 <- coord_shift(df_coord_count_robinson, 100000)
df_shift3 <- coord_shift(df_coord_count_robinson, 190000)
dotsize <- 1
trisize <- 0.4


# 20240907. Figure. Data collection map - Final (hopefully)
plot_europe <- ggplot(data = world_robinson) +
  geom_sf(color='black', aes(fill=log10(n)), show.legend = F) + # Contamination
  geom_sf(data=subset(world_robinson, is.na(n)),fill='#EEEEEE', color='black') + # Empty points
  geom_sf(data=dispute_robinson, fill='#BBBBBB', color=NA, linetype=2) + # Border dispute
  geom_sf(data=subset(df_shift1, n<100), fill='black', shape=25, size=trisize) + # Institute locations
  geom_sf(data=subset(df_shift2, n<100), color='black', aes(fill=log10(n)), shape=21, size=dotsize, stroke=0.5, show.legend = F) + # Institute locations
  geom_sf(data=subset(df_shift1, n>=100), fill='black', shape=25, size=trisize) + # Institute locations
  geom_sf(data=subset(df_shift2, n>=100), color='black', aes(fill=log10(n)), shape=21, size=dotsize, stroke=0.5, show.legend = F) + # Institute locations
  scale_fill_gradient2(low='#BBBBBB',mid='#F3CEBB',high='#B31312', midpoint=2.2) +
  xlim(c(-1000000, 3000000)) +
  ylim(c(3600000, 7400000)) +
  theme_void() +
  theme(panel.background = element_rect(fill='white', color='black'))

ggplot(data = world_robinson) +
  geom_sf(color='black', aes(fill=log10(n))) + # Contamination
  geom_sf(data=subset(world_robinson, is.na(n)),fill='#EEEEEE', color='black') + # Empty points
  geom_sf(data=dispute_robinson, fill='#BBBBBB', color=NA, linetype=2) + # Border dispute
  geom_sf(data=subset(df_shift2, n<100), fill='black', shape=25, size=trisize) + # Institute locations
  geom_sf(data=subset(df_shift3, n<100), color='black', aes(fill=log10(n)), shape=21, size=dotsize, stroke=0.5) + # Institute locations
  geom_sf(data=subset(df_shift2, n>=100), fill='black', shape=25, size=trisize) + # Institute locations
  geom_sf(data=subset(df_shift3, n>=100), color='black', aes(fill=log10(n)), shape=21, size=dotsize, stroke=0.5) + # Institute locations
  scale_fill_gradient2(low='#BBBBBB',mid='#F3CEBB',high='#B31312', midpoint=2.2, name='log10(Samples)') +
  annotation_custom(
    grob = ggplotGrob(plot_europe),
    xmin = -4000000,
    xmax = +5000000,
    ymin = -7000000,
    ymax = +3000000
  ) +
  annotate("rect",fill=NA, color='black', linewidth=0.5,
           xmin=-1000000, xmax=3000000,
           ymin=3600000, ymax=7500000)+
  annotate("rect",fill=NA, color='black', linewidth=0.5,
           xmin=-4000000, xmax=5000000,
           ymin=-6270000, ymax=2260000)+
  annotate("segment", colour = "black", linetype=2,
           x = -4000000, xend = -1000000, 
           y = 2260000, yend = 3600000) + # Left
  annotate("segment", colour = "black", linetype=2,
           x = 3000000, xend = 5000000, 
           y = 3600000, yend = 2260000) + # Right
  theme_void() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.95,0.7))

ggsave(filename=paste0(path_root, '/pdfs/20241024_FigureX_data_collection_map.pdf'), unit="in",width=10.0, height=4.0)

