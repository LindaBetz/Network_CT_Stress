
# ---------------------------- 1: Load libraries & packages -----------------------------
library(qgraph)
library(igraph)
library(bootnet)
library(coin)
library(purrr)
library(dplyr)

# all data sets are available at https://www.icpsr.umich.edu/icpsrweb/

# original sample
biomarker_data_original <- da29282.0001 # CTQ, PSS in here
demographic_data_original <-
  da04652.0001 # demographic & clinical vars in here

# refresher sample
biomarker_data_refresher <- da36901.0001 # CTQ, PSS in here
demographic_data_refresher <-
  da36532.0001 # demographic & clinical vars in here

# ---------------------------- 2: Data wrangling & sample characteristics -----------------------------

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
  # pos
  "Things were not going your way",
  # pos
  "Unable to control irritations in life",
  # pos
  "Did not feel on top of things"
  # pos
)

# filter those people with no missing values
relevant_IDs_original <- biomarker_data_original %>%
  select(.,
         matches("M2ID|B4QCT_EA|B4QCT_SA|B4QCT_PA|B4QCT_EN|B4QCT_PN|B4Q4")) %>%
  mutate(na_per_row = rowSums(is.na(.) / 15)) %>% #M2ID never missing
  filter(na_per_row <= 13 / 15) %>%
  transmute(M2ID)


nrow(relevant_IDs_original) # 1252 ==> 3 people do not have at least two variables available, we exclude them


# create data set to calculate sample statstiscs
desc_data_original <- biomarker_data_original %>%
  left_join(demographic_data_original, by = "M2ID") %>%
  filter(M2ID %in% relevant_IDs_original$M2ID) %>%
  transmute(
    Age = B1PAGE_M2.y,
    Sex = ifelse(B1PRSEX.y == "(1) MALE", 1, 0),
    Site =  B4ZSITE,
    MIDUS_Sample = SAMPLMAJ.x,
    CESD = B4QCESD,
    PSS = B4QPS_PS,
    Emotional_Abuse = B4QCT_EA,
    Sexual_Abuse = B4QCT_SA,
    Physical_Abuse = B4QCT_PA,
    Emotional_Neglect = B4QCT_EN,
    Phyiscal_Neglect = B4QCT_PN,
    
    Ethnicity = ifelse(
      B1PF7A == "(1) WHITE",
      "White",
      ifelse(
        B1PF7A == "(2) BLACK AND/OR AFRICAN AMERICAN",
        "African-American",
        "Other"
      )
    ),
    Education = ifelse(
      as.numeric(B1PB1) <= 3,
      "less than high school",
      ifelse(
        as.numeric(B1PB1) > 3 &
          as.numeric(B1PB1) <= 8,
        "graduated at least high school or obtained GED",
        ifelse(as.numeric(B1PB1) > 8, "4-year college degree or more", NA)
      )
    )
  ) %>% mutate(Sample = "Original")

# frequency/count data
desc_data_original %>%
  select(., c(Site, MIDUS_Sample, Ethnicity, Education)) %>%
  map(~ table(.) / sum(!is.na(.))) %>%
  map(~ round(., 3))

# continous data
desc_data_original %>%
  select(
    .,
    c(
      Age,
      Sex,
      CESD,
      PSS,
      Emotional_Abuse,
      Sexual_Abuse,
      Physical_Abuse,
      Emotional_Neglect,
      Phyiscal_Neglect
    )
  ) %>%
  summarise_all(c("mean", "sd"), na.rm = TRUE) %>%
  round(., 1)

# make a data set to be used in the estimation of the network
graph_data_original <- biomarker_data_original %>%
  # filter(SAMPLMAJ == "(01) MAIN RDD") %>%
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
# ---------------------------- 3: Network estimation & visualization -----------------------------

# estimate network
graph_original <- estimateNetwork(
  graph_data_original,
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

# estimate communities via walktrap
g <- as.igraph(qgraph(graph_original$graph), attributes = TRUE)
wtc <- walktrap.community(g)
wtc$membership


# plot network
graph_original_plot <- qgraph(
  graph_original$graph,
  layout = "spring",
  theme = "Borkulo",
  layout.par = list(max.delta = 10),
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  GLratio = 1.9,
  groups =  recode(
    wtc$membership,
    `2` = "Childhood Trauma",
    `1` = "Perceived Helplessness",
    `3` = "Perceived Self-Efficacy"
  ),
  layoutOffset = c(-0.26, 0),
  layoutScale = c(0.9, 1),
  label.cex = 1.25,
  filetype = "png",
  color =  recode(
    wtc$membership,
    `2` = "#D3D3D3",
    `1` = "#FF9770",
    `3` = "#5DBCD2"
  ),
  nodeNames = colnames(graph_original$graph),
  filename = "main_figure",
  legend.cex = 0.69
)


# ---------------------------- 4: Sensitvity Analyses -----------------------------
# case-drop bootstrapping
results_case <-
  bootnet(graph_original, type = "case", nBoots = 1000)
corStability(results_case) # 0.75 for edge & strength

# supplementary plot: case-dropping strength
#png(filename="case_dropping_strength.png", width=600, height=400)
plot(results_case, statistics = "strength") +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
  theme(
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
#dev.off()

# supplementary plot: case-dropping edge
#png(filename="case_dropping_edge.png", width=600, height=400)
plot(results_case, statistics = c("edge")) +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1), limits = c(0.5, 1)) +
  theme(
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
#dev.off()


# regular bootstrapping for edge weights
results_boot <- bootnet(graph_original, nBoots = 1000)

# supplementary plot: bootstrapped edges
#png("bootstrapped_edges.png", height=1400, width=600, pointsize = 12.5)
plot(
  results_boot,
  order = "sample",
  split0 = TRUE,
  labels = TRUE,
  legend = TRUE,
  prop0_cex = 2
) + theme(text = element_text(size = 13))
#dev.off()

# ---------------------------- 5: Replication Analyses -----------------------------
# filter those people with no missing values
relevant_IDs_refresher <- biomarker_data_refresher %>%
  select(.,
         matches(
           "MRID|RA4QCT_EA|RA4QCT_SA|RA4QCT_PA|RA4QCT_EN|RA4QCT_PN|RA4Q4"
         )) %>%
  mutate(na_per_row = rowSums(is.na(.) / 15)) %>% #MRID never missing
  filter(na_per_row <= 13 / 15) %>%
  transmute(MRID)

nrow(relevant_IDs_refresher) # 862 ==> one person does not have at least two variables, we thus exclude

# create data set to calculate sample statstiscs
desc_data_refresher <- biomarker_data_refresher %>%
  left_join(baseline_refresher_data, by = "MRID") %>%
  filter(MRID %in% relevant_IDs_refresher$MRID) %>%
  transmute(
    Age = RA1PRAGE.y,
    Sex = ifelse(RA1PRSEX.y == "(1) MALE", 1, 0),
    Site =  RA4ZSITE,
    MIDUS_Sample = SAMPLMAJ.y,
    CESD = RA4QCESD,
    PSS = RA4QPS_PS,
    Emotional_Abuse = RA4QCT_EA,
    Sexual_Abuse = RA4QCT_SA,
    Physical_Abuse = RA4QCT_PA,
    Emotional_Neglect = RA4QCT_EN,
    Phyiscal_Neglect = RA4QCT_PN,
    
    Ethnicity = ifelse(
      RA1PF7A == "(1) WHITE",
      "White",
      ifelse(
        RA1PF7A == "(2) BLACK AND/OR AFRICAN AMERICAN",
        "African-American",
        "Other"
      )
    ),
    Education = ifelse(
      as.numeric(RA1PB1) <= 3,
      "less than high school",
      ifelse(
        as.numeric(RA1PB1) > 3 &
          as.numeric(RA1PB1) <= 8,
        "graduated at least high school or obtained GED",
        ifelse(as.numeric(RA1PB1) > 8, "4-year college degree or more", NA)
      )
    )
  ) %>% mutate(Sample = "Refresher")


# frequency/count data
desc_data_refresher %>%
  select(., c(Ethnicity, Education)) %>%
  map(~ table(.) / sum(!is.na(.))) %>%
  map(~ round(., 3))

# continous data
desc_data_refresher %>%
  select(
    .,
    c(
      Age,
      Sex,
      CESD,
      PSS,
      Emotional_Abuse,
      Sexual_Abuse,
      Physical_Abuse,
      Emotional_Neglect,
      Phyiscal_Neglect
    )
  ) %>%
  summarise_all(c("mean", "sd"), na.rm = TRUE) %>%
  round(., 1)

# statistical comparison of samples

combined_data <-
  rbind.data.frame(desc_data_original, desc_data_refresher)

combined_data %>%
  select(., c(Sex, Ethnicity, Education)) %>%
  map(
    ~ chisq_test(
      as.factor(.) ~ as.factor(combined_data$Sample),
      data = combined_data,
      distribution = "approximate"
    )
  )

combined_data %>%
  select(
    .,
    c(
      Age,
      CESD,
      PSS,
      Emotional_Abuse,
      Sexual_Abuse,
      Physical_Abuse,
      Emotional_Neglect,
      Phyiscal_Neglect
    )
  ) %>%
  map(
    ~ oneway_test(
      as.numeric(.) ~ as.factor(combined_data$Sample),
      data = combined_data,
      distribution = "approximate"
    )
  )



# extract relevant variables from data set, basic "preprocessing" as above
graph_data_refresher <- biomarker_data_refresher %>%
  # filter(SAMPLMAJ == "(01) MAIN RDD") %>%
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

# estimate refresher network
graph_refresher <- estimateNetwork(
  graph_data_refresher,
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

# basic comparison: correlate the two partial correlation matrices
cor(c(graph_original$graph), c(graph_refresher$graph)) # 0.9242128

# Network Comparison for Differences in Structure, Global Strength & Individual Edges
set.seed(1337)
# may take longer to run on a standard PC!
compare_12 <-
  NCT(
    graph_original,
    graph_refresher,
    it = 1000,
    test.edges = TRUE,
    edges = 'all',
    progressbar = TRUE
  )

# NCT Structure differences
compare_12$nwinv.pval #  0.783

# NCT quantification of differences: count significantly different edges (total number of edges 105)
sum(compare_12$einv.pvals$"p-value" < 0.05) # 0
sum(p.adjust(compare_12$einv.pvals$"p-value", "BH") < 0.05) # 0

# NCT global strength
compare_12$glstrinv.pval # 0.496

# correlate strength centrality indices
graph_original_strength <-
  centralityTable(graph_original$graph) %>% filter(measure == "Strength") %>%
  transmute(value)

graph_refresher_strength <-
  centralityTable(graph_refresher$graph) %>% filter(measure == "Strength") %>%
  transmute(value)

cor(graph_original_strength, graph_refresher_strength) # 0.9562717

# plot refresher network & combined network
# estimate communities via walktrap
walktrap.community(as.igraph(qgraph(graph_refresher$graph), attributes = TRUE))$membership

# estimate combined network
graph_combined <- estimateNetwork(
  rbind.data.frame(graph_data_original, graph_data_refresher),
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

# walktrap combined
walktrap.community(as.igraph(qgraph(graph_combined$graph), attributes = TRUE))$membership

# communities are the same across all three graphs. That's why we use the original object for grouping in the plots


png(width=1200, height=450, "combined_plot.png")
par(mfrow=c(1,3))
qgraph(
  graph_original$graph,
  layout = graph_original_plot$layout,
  theme = "Borkulo",
  labels = TRUE,
  layout.par = list(max.delta = 10),
  labels = c(1:15),
  pie = rep(1, 15),
  title="a) Original",
  title.cex=3,
  vsize = 14,
  pieBorder = 0.3,
  pieColor = recode(
    wtc$membership,
    `2` = "#FFD670",
    `1` = "#FF9770",
    `3` = "#D3D3D3",
    
  ),
  groups = c(rep("1) Childhood Trauma", 5), rep("2) Perceived Stress", 10)),
  label.cex = 2,
  color = c(rep("#D3D3D3", 5),
            rep("#5DBCD2", 10)),
  nodeNames = colnames(graph_data),
  legend=FALSE,
  maximum =0.4843399,
  minimum = 0
  
)

qgraph(
  graph_refresher$graph,
  layout = graph_original_plot$layout,
  theme = "Borkulo",
  labels = TRUE,
  layout.par = list(max.delta = 10),
  labels = c(1:15),
  pie = rep(1, 15),
  title="b) Replication",
  title.cex=3,
  vsize = 14,
  pieBorder = 0.3,
  pieColor = recode(
    wtc$membership,
    `2` = "#FFD670",
    `1` = "#FF9770",
    `3` = "#D3D3D3",
    
  ),
  groups = c(rep("1) Childhood Trauma", 5), rep("2) Perceived Stress", 10)),
  label.cex = 2,
  color = c(rep("#D3D3D3", 5),
            rep("#5DBCD2", 10)),
  nodeNames = colnames(graph_data),
  legend=FALSE,
  maximum = 0.4843399,
  minimum = 0
  

)

qgraph(
  graph_combined$graph,
  layout = graph_original_plot$layout,
  theme = "Borkulo",
  labels = TRUE,
  layout.par = list(max.delta = 10),
  labels = c(1:15),
  title="c) Combined",
  title.cex=3,
  vsize = 14,
  pie = rep(1, 15),
  pieBorder = 0.3,
  pieColor = recode(
    wtc$membership,
    `2` = "#FFD670",
    `1` = "#FF9770",
    `3` = "#D3D3D3",
    
  ),
  groups = c(rep("1) Childhood Trauma", 5), rep("2) Perceived Stress", 10)),
  label.cex = 2,
  color = c(rep("#D3D3D3", 5),
            rep("#5DBCD2", 10)),
  nodeNames = colnames(graph_data),
  legend=FALSE,
  maximum = 0.4843399,
  minimum = 0
  
  
)
dev.off()
