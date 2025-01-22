# institute_network
# Responsible for...
#   Institute network

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
library(igraph)
library(ggraph)
library(ggforce)
library(concaveman)

# Read data
path_root <- "../../"
path_institute_agg <- paste0(path_root, '/analysis/4_check_misidentified/misidentified_institute.txt') 

df_institute_agg <- fread(path_institute_agg, header=F) %>% as.data.frame()
names(df_institute_agg) <- c("RorID", "Institute", "Country", "_3", "_4", "_5", "Frac_institute",
                             "match", "mismatch", "total", "Frac_aggregation", 'Parents')

df_institute_agg_filt <- subset(df_institute_agg, total >= 10 & RorID != "DISCARD")
nrow(df_institute_agg_filt) # 1,179 target institutes


# Build graph
df_graph <- df_institute_agg_filt[,c('RorID','Parents', 'total', 'Frac_aggregation')]
df_graph$label <- ifelse(df_graph$RorID == df_graph$Parents, df_graph$RorID, NA)
df_graph$color <- cut(df_graph$Frac_aggregation,
                      breaks = c(0.00, 0.000000000001, 0.05, 0.10, 0.20, 1.00),
                      labels = c('#337357', '#B3B800', '#E59500','#B22222', '#5A0000'),
                      include.lowest = T)
  
my_nodes <- df_graph[,c('RorID', 'total', 'color', 'label', 'Frac_aggregation')]
my_edges <- df_graph[,c('RorID','Parents')] %>% 
  separate_rows(Parents, sep='; ') %>%
  subset(RorID != Parents)

g <- graph_from_data_frame(my_edges, directed=T, vertices = my_nodes)


# Assign group by progenitor (which is country)
root_nodes <- V(g)[!is.na(V(g)$label)]
find_descendants <- function(graph, root) {
  subgraph <- subcomponent(graph, root, mode = "in")  # Find all descendants
  return(subgraph)
}

root_clusters <- list()
for (root in root_nodes) {
  descendants <- find_descendants(g, root)
  root_clusters[[length(root_clusters) + 1]] <- descendants
}

V(g)$group <- NA  # Initialize a group vector with NA

for (i in seq_along(root_clusters)) {
  for (v in root_clusters[[i]]) {
    # Assign only if the node doesn't already have a group
    # There are institutes with multiple affiliation, we choose one randomly
    if (is.na(V(g)$group[v])) {  
      V(g)$group[v] <- i
    }
  }
}

V(g)$group <- as.factor(V(g)$group)
V(g)$color <- factor(V(g)$color, levels = c('#337357', '#B3B800', '#E59500','#B22222', '#5A0000'))


# setup weights
weights <- rep(0.20, length(E(g)))
for (edge in E(g)) {
  group <- V(g)[ends(g, edge)[2]]$group
  if (group %in% c("2", "3")) { # China and USA
    weights[edge] <- 0.65
  }
  # weights[edge] <- group_size[group,'count'] / max(group_size[,'count'])
}

# Make initial layout
set.seed(15) # was 12.. now changed
# 10,14,15 is decent

layout <-  layout_with_fr(g, niter = 3000, grid='nogrid',
                          weights = weights)


# Find subgraph where size is 2. this is because geom_mark_hull doesn't mark group size of 2
group_attr <- V(g)$group
group_size <- aggregate(x = list(count = group_attr), by = list(group = group_attr), FUN = length)
groups_of_size_2 <- group_size[group_size$count <= 2, "group"]
vertices_of_interest <- V(g)[group_attr %in% groups_of_size_2]
subgraph <- induced_subgraph(g, vids = vertices_of_interest)
subgraph_layout <- layout[as.numeric(vertices_of_interest), ] %>% 
  as.data.frame() %>% 
  cbind(vertices_of_interest$group)
names(subgraph_layout) <- c('x','y', 'group')


# Figure. institute network. sample>=10.
# Color is intended to be green2red, but fixed to grey2red. won't fix this bad code
palette <- c('#CCCCCC', '#FFD23F', '#FFA500', '#D73027', '#5A0000')
colormapper <- c('#337357'= palette[1],
                 '#B3B800'= palette[2],
                 '#E59500'= palette[3],
                 '#B22222'= palette[4],
                 '#5A0000'= palette[5]) 


ggraph(g, layout=layout) +
  geom_edge_link(color='#777777') +
  geom_node_point(aes(size=log10(total)*20, fill=color), 
                   color='black', shape=21, show.legend = T) +
  geom_node_label(aes(label = label), color=NA, repel=F,
                  vjust = 0.2, hjust=0.05, size = 3, fill='white', alpha=0.5) +
  geom_node_text(aes(label = label), color='black', repel=F, fontface=1,
                 vjust = 0.0, hjust=0, size = 3, alpha=1) +
  geom_mark_hull(aes(x=x, y=y, group=group),
                 expand=unit(3, "mm"), radius = unit(3, "mm"), concavity = 3,) +
  geom_mark_hull(data=subgraph_layout,
                 aes(x=x, y=y, group=group),
                 expand=unit(2, "mm"), radius = unit(3, "mm"), concavity = 3,) +
  scale_fill_manual(name="Contaminated (%)",
                    values=colormapper,
                    labels=c('#337357'="=0%", '#B3B800'=">0%",
                             '#E59500'=">5%", '#B22222'=">10%", 
                             '#5A0000'=">20%")) +
  scale_size_continuous(name = "Data size",
                        breaks = log10(c(10, 100, 1000, 10000)) * 20,
                        labels = c("10", "100", "1,000", "10,000")) +
  theme_void()

ggsave(filename=paste0(path_root, '/pdfs/202411111_FigureX_institute_network.pdf'), unit="in",width=9.5, height=8.0)


# Count terminal nodes
v_terminal <- V(g_tmp)[degree(g_tmp, mode = "in") == 0]

length(v_terminal) # total 881
sum(v_terminal$Frac_aggregation != 0.00) # > 3%: 285
sum(v_terminal$Frac_aggregation > 0.03) # > 3%: 285
sum(v_terminal$Frac_aggregation > 0.20) # > 20%: 141

table(v_terminal$color) # 0%: 596, 3%:20, 5%:64, 10%: 60, 20%: 141
