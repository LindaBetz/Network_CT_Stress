# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Code for    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                                               #
#         Relationships between childhood trauma and subjective experiences of stress           #
#                     in the general population: a network perspective.                         #
#                                 developed by L. Betz                                          #
#                                                                                               #
#                           - Analysis reported in supplementary material -                     #
#                                                                                               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# ---------------------------------- 1: Load libraries & packages ----------------------------------
library(qgraph)
library(igraph)
library(bootnet)
library(NetworkComparisonTest)
library(EGAnet)
library(purrr)
library(dplyr)

# all data sets are available at https://www.icpsr.umich.edu/icpsrweb/

# original sample (= Biomarker "original")
biomarker_data_original <- da29282.0001 # CTQ, PSS in here
demographic_data_original <-
  da04652.0001 # demographic & clinical vars in here

# replication sample (= Biomarker refresher)
biomarker_data_replication <- da36901.0001 # CTQ, PSS in here
demographic_data_replication <-
  da36532.0001 # demographic & clinical vars in here

# ------------------------------------- 2: Data preparation --------------------------------------
# variable names. Note that for the PSS, we reworded the positive variables for visualization to make interpretation easier
var_names <- c(
  "Upset by something unexpected",
  "Unable to control important things",
  "Felt nervous and stressed",
  "Not confident to handle personal problems",
  # pos
  "Things were not going your way",
  # pos
  "Could not cope with all things to do",
  "Unable to control irritations in life",
  # pos
  "Did not feel on top of things",
  # pos
  "Angered by things outside control",
  "Difficulties piling up can't overcome",
  "Emotional Abuse",
  "Physical Abuse",
  "Sexual Abuse",
  "Emotional Neglect",
  "Physical Neglect"
)

# names of PSS-variables that will be recoded
recode_vars <- c(
  "Not confident to handle personal problems",
  "Things were not going your way",
  "Unable to control irritations in life",
  "Did not feel on top of things"
  
)

## .......................... Original Sample  ..........................

### filter those people with no missing values
relevant_IDs_original <- biomarker_data_original %>%
  select(.,
         matches("M2ID|B4QCT_EA|B4QCT_SA|B4QCT_PA|B4QCT_EN|B4QCT_PN|B4Q4")) %>%
  mutate(na_per_row = rowSums(is.na(.) / 15)) %>% #M2ID never missing
  filter(na_per_row <= 13 / 15) %>% # at least two variables available
  transmute(M2ID)

nrow(relevant_IDs_original) # 1252 ==> 3 people do not have at least two variables available, we exclude them

### make a data set to be used in the estimation of the network
graph_data_original <- biomarker_data_original %>%
  select(.,
         matches("B4QCT_EA|B4QCT_SA|B4QCT_PA|B4QCT_EN|B4QCT_PN|B4Q4")) %>%
  `colnames<-`(var_names) %>%
  mutate_all(as.numeric) %>%
  mutate_at(recode_vars,
            ~ recode(
              # recode positive items
              .,
              `1` = 5,
              `2` = 4,
              `3` = 3,
              `4` = 2,
              `5` = 1,
              .missing = NA_real_
            )) %>%
  select(
    # change order of items, to make plot nicer later
    `Emotional Neglect`,
    `Physical Neglect`,
    `Emotional Abuse`,
    `Physical Abuse`,
    `Sexual Abuse`,
    everything()
  )


graph_original <- estimateNetwork(
  graph_data_original,
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

## ..........................  Replication Sample  ..........................

# extract relevant variables from data set, basic "preprocessing" as above
graph_data_replication <- biomarker_data_replication %>%
  select(.,
         matches("RA4QCT_EA|RA4QCT_SA|RA4QCT_PA|RA4QCT_EN|RA4QCT_PN|RA4Q4")) %>%
  `colnames<-`(var_names) %>%
  mutate_all(as.numeric) %>%
  mutate_at(recode_vars,
            ~ recode(
              # recode positive items
              .,
              `1` = 5,
              `2` = 4,
              `3` = 3,
              `4` = 2,
              `5` = 1,
              .missing = NA_real_
            )) %>%
  select(
    # change order of items, to make plot nicer later
    `Emotional Neglect`,
    `Physical Neglect`,
    `Emotional Abuse`,
    `Physical Abuse`,
    `Sexual Abuse`,
    everything()
  )

graph_replication <- estimateNetwork(
  graph_data_replication,
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

# ------------------------------------- 3: Robustness Analyses --------------------------------------

## .......................... original sample ..........................

### _____________ case-drop bootstrapping _____________
results_case_original <-
  bootnet(graph_original, type = "case", nBoots = 1000)
corStability(results_case) # 0.75 for edge & strength

#### supplementary plot: case-dropping strength
tiff(filename = "case_dropping_strength_original.tiff",
     width = 600,
     height = 400)
plot(results_case_original, statistics = "strength") +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
  theme(
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
dev.off()

#### supplementary plot: case-dropping edge
tiff(filename = "case_dropping_edge_original.tiff",
     width = 600,
     height = 400)
plot(results_case_original, statistics = c("edge")) +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
  theme(
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
dev.off()


### _____________ regular bootstrapping for edge weights _____________
results_boot_original <- bootnet(graph_original, nBoots = 1000)

#### supplementary plot: bootstrapped edges
tiff(
  "bootstrapped_edges_original.tiff",
  height = 1400,
  width = 600,
  pointsize = 13
)
plot(
  results_boot_original,
  order = "sample",
  split0 = TRUE,
  labels = TRUE,
  legend = TRUE,
  prop0_cex = 2
) + theme(text = element_text(size = 13))
dev.off()

## .......................... replication sample ..........................
### _____________ case-drop bootstrapping _____________
results_case_replication <-
  bootnet(graph_replication, type = "case", nBoots = 1000)
corStability(results_case_replication) # 0.75 for edge & strength

#### supplementary plot: case-dropping strength
tiff(filename = "case_dropping_strength_replication.tiff",
     width = 600,
     height = 400)
plot(results_case_replication, statistics = "strength") +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
  theme(
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
dev.off()

#### supplementary plot: case-dropping edge
tiff(filename = "case_dropping_edge_replication.tiff",
     width = 600,
     height = 400)
plot(results_case_replication, statistics = c("edge")) +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
  theme(
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
dev.off()


### _____________ regular bootstrapping for edge weights _____________
results_boot_replication <-
  bootnet(graph_replication, nBoots = 1000)

# supplementary plot: bootstrapped edges
tiff(
  "bootstrapped_edges_replication.tiff",
  height = 1400,
  width = 600,
  pointsize = 13
)
plot(
  results_boot_replication,
  order = "sample",
  split0 = TRUE,
  labels = TRUE,
  legend = TRUE,
  prop0_cex = 2
) + theme(text = element_text(size = 13))
dev.off()


# ------------------------------------- 4: visualization of original, replication network & combined network --------------------------------------
### _____________ estimate communities via walktrap for original sample _____________
wtc <-
  walktrap.community(as.igraph(qgraph(graph_original$graph), attributes = TRUE))



### _____________ estimate communities via walktrap for replication sample _____________
walktrap.community(as.igraph(qgraph(graph_replication$graph), attributes = TRUE))$membership

### _____________ estimate combined network _____________
graph_combined <- estimateNetwork(
  rbind.data.frame(graph_data_original, graph_data_replication),
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

### _____________ estimate walktrap for combined sample  _____________
walktrap.community(as.igraph(qgraph(graph_combined$graph), attributes = TRUE))$membership

# ==> communities are the same across all three graphs. That's why we use the original object for grouping in the plots


### _____________ layout for network  _____________
layout_network <- as.matrix(data.frame(
  x =  c(
    0.702164557,
    0.980149473,
    0.202511180,
    0.579496614,
    0.500319086,
    0.382166936,
    0.265510012,
    0.000000000,
    1.000000000,
    0.717622685,
    0.270821989,
    0.917685999,
    0.715207072,
    0.004062325,
    0.440950384
  ),
  y = c(
    0.88338016,
    0.80279178,
    0.78809824,
    0.63657324,
    1.00000000,
    0.09806150,
    0.33345000,
    0.37831538,
    0.32678846,
    0.05273452,
    0.61391449,
    0.00000000,
    0.41124686,
    0.04156892,
    0.24338861
  )
))


### _____________ supplementary plot: combined, original and replication network next to each other
tiff(width = 1200, height = 450, "combined_plot.tiff")
par(mfrow = c(1, 3))
qgraph(
  graph_original$graph,
  layout = layout_network,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  legend = F,
  GLratio = 1.1,
  groups =  recode(
    wtc$membership,
    `2` = "Childhood Trauma",
    `1` = "Perceived Helplessness",
    `3` = "Perceived Self-Efficacy"
  ),
  layoutOffset = c(-0.05, 0),
  layoutScale = c(1, 1.05),
  label.cex = 0.99,
  color =  c("grey",
             "#EBCC2A",
             "#78B7C5"),
  label.prop = 0.96,
  vsize = 13,
  DoNotPlot = F,
  nodeNames = colnames(graph_original$graph),
  minimum = 0,
  maximum = 0.483,
  title = "a) Original",
  title.cex = 3
)

qgraph(
  graph_replication$graph,
  layout = layout_network,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  legend = F,
  GLratio = 1.1,
  groups =  recode(
    wtc$membership,
    `2` = "Childhood Trauma",
    `1` = "Perceived Helplessness",
    `3` = "Perceived Self-Efficacy"
  ),
  layoutOffset = c(-0.05, 0),
  layoutScale = c(1, 1.05),
  label.cex = 0.99,
  color =  c("grey",
             "#EBCC2A",
             "#78B7C5"),
  label.prop = 0.96,
  vsize = 13,
  DoNotPlot = F,
  nodeNames = colnames(graph_original$graph),
  minimum = 0,
  maximum = 0.483,
  title = "b) Replication",
  title.cex = 3
)

qgraph(
  graph_combined$graph,
  layout = layout_network,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  legend = F,
  GLratio = 1.1,
  groups =  recode(
    wtc$membership,
    `2` = "Childhood Trauma",
    `1` = "Perceived Helplessness",
    `3` = "Perceived Self-Efficacy"
  ),
  layoutOffset = c(-0.05, 0),
  layoutScale = c(1, 1.05),
  label.cex = 0.99,
  color =  c("grey",
             "#EBCC2A",
             "#78B7C5"),
  label.prop = 0.96,
  vsize = 14,
  DoNotPlot = F,
  nodeNames = colnames(graph_original$graph),
  minimum = 0,
  maximum = 0.483,
  title = "c) Combined",
  title.cex = 3
  
)

dev.off()

# ------------------------------------- 5: Alternative community detection  --------------------------------------

g <-
  as.igraph(qgraph(graph_original$graph, DoNotPlot = TRUE), attributes = TRUE)

wtc <- walktrap.community(g)
set.seed(234)
sgc <- spinglass.community(g)

wtc$membership == sgc$membership # results are identical (walktrap vs. spinglass)

ega <- EGA(graph_data_original, plot.EGA = F)

recode(ega$wc,
       `1` = 2,
       `2` = 1,
       `3` = 3) == wtc$membership # results are identical (walktrap vs. EGA)


# ------------------------------------- 6: Centrality plots  --------------------------------------
centralityPlot(qgraph(graph_original$graph), orderBy = "Strength")
centralityPlot(qgraph(graph_replication$graph), orderBy = "Strength")


