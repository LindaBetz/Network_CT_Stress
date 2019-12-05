
# ---------------------------- 1: Load libraries & packages -----------------------------
library(qgraph)
library(igraph)
library(bootnet)
library(NetworkComparisonTest)
library(coin)
library(purrr)
library(dplyr)

# all data sets are available at https://www.icpsr.umich.edu/icpsrweb/

# original sample (= Biomarker "original")
biomarker_data_original <- da29282.0001 # CTQ, PSS in here
demographic_data_original <- da04652.0001 # demographic & clinical vars in here

# replication sample (= Biomarker refresher)
biomarker_data_replication <- da36901.0001 # CTQ, PSS in here
demographic_data_replication <- da36532.0001 # demographic & clinical vars in here

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

# ----------------------------------- ORIGINAL SAMPLE 

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
    Age = B1PAGE_M2.x,
    Sex = ifelse(B1PRSEX.x == "(1) MALE", 1, 0),
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
         matches("B1PRSEX|B4QCT_EA|B4QCT_SA|B4QCT_PA|B4QCT_EN|B4QCT_PN|B4Q4")) %>%
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
    Sex = ifelse(B1PRSEX == "(1) MALE", 1, 0),
    everything()
  )

# ------------------------------------ REPLICATION SAMPLE

# filter those people with no missing values
relevant_IDs_replication <- biomarker_data_replication %>%
  select(.,
         matches(
           "MRID|RA4QCT_EA|RA4QCT_SA|RA4QCT_PA|RA4QCT_EN|RA4QCT_PN|RA4Q4"
         )) %>%
  mutate(na_per_row = rowSums(is.na(.) / 15)) %>% #MRID never missing
  filter(na_per_row <= 13 / 15) %>%
  transmute(MRID)

nrow(relevant_IDs_replication) # 862 ==> one person does not have at least two variables, we thus exclude


# create data set to calculate sample statstiscs
desc_data_replication <- biomarker_data_replication %>%
  left_join(demographic_data_replication, by = "MRID") %>%
  filter(MRID %in% relevant_IDs_replication$MRID) %>%
  transmute(
    Age = RA1PRAGE.x,
    Sex = ifelse(RA1PRSEX.x == "(1) MALE", 1, 0),
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
  ) %>% mutate(Sample = "Replication")


# frequency/count data
desc_data_replication %>%
  select(., c(Ethnicity, Education)) %>%
  map(~ table(.) / sum(!is.na(.))) %>%
  map(~ round(., 3))

# continous data
desc_data_replication %>%
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

# ----------------- Statistical Comparison
combined_data <-
  rbind.data.frame(desc_data_original, desc_data_replication)

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

# ------------------------------------ COMBINED SAMPLE
# descriptives for combined sample

desc_combined_data <-
  rbind.data.frame(desc_data_original, desc_data_replication)

# frequency/count data
desc_combined_data %>%
  select(., c(Ethnicity, Education)) %>%
  map(~ table(.) / sum(!is.na(.))) %>%
  map(~ round(., 3))

# continous data
desc_combined_data %>%
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

# ------------------------------------ MALES vs. FEMALES
# ----------------- Descriptives

# 0 = women, 1 = men
# frequency/count data
desc_combined_data %>%
  select(., c(Sex, Ethnicity, Education)) %>%
  split(.$Sex) %>% 
  map(~select(., -c("Sex"))) %>% # remove sex variable after grouping
  map(~ c(table(.$Ethnicity) / sum(!is.na(.$Ethnicity)),
         table(.$Education) / sum(!is.na(.$Education)))) %>%
  map(~ round(., 3))

# continous data
desc_combined_data  %>%
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
  ) %>%  group_by(Sex) %>% # males
  summarise_all(c("mean", "sd"), na.rm = TRUE) %>%
  round(., 1)

# ----------------- Statistical Comparison
set.seed(1)
desc_combined_data %>%
  select(., c(Ethnicity, Education)) %>%
  map(
    ~ chisq_test(
      as.factor(.) ~ as.factor(desc_combined_data$Sex),
      data = combined_data,
      distribution = "approximate"
    )
  )

desc_combined_data %>%
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
      as.numeric(.) ~ as.factor(desc_combined_data$Sex),
      data = combined_data,
      distribution = "approximate"
    )
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
g <- as.igraph(qgraph(graph_original$graph, DoNotPlot=TRUE), attributes = TRUE)
wtc <- walktrap.community(g)

# layout for network
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

# plot network
graph_original_plot <- qgraph(
  graph_original$graph,
  layout = layout_network,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  GLratio = 1.1,
  groups =  recode(
    wtc$membership,
    `2` = "Childhood Trauma",
    `1` = "Perceived Helplessness",
    `3` = "Perceived Self-Efficacy"
  ),
  layoutOffset = c(-0.05, 0),
  layoutScale = c(1.14, 1.05),
  label.cex = 0.99,
  filetype = "png",
  legend.mode="style1",
  color =  c("grey",
    "#EBCC2A",
    "#78B7C5"
  ),
  label.prop = 0.96,
  vsize=5.1,
  DoNotPlot=F,
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

# ---------------------------- 5: Comparison Sex Differences -----------------------------

# here, we first merge the original and replication sample to retain sufficient power
graph_data_sex <- biomarker_data_original %>%
  filter(M2ID %in% relevant_IDs_original$M2ID) %>%
  select(.,
         matches("B1PRSEX|B4QCT_EA|B4QCT_SA|B4QCT_PA|B4QCT_EN|B4QCT_PN|B4Q4")) %>%
  `colnames<-`(c("Sex", var_names)) %>%
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
  select(Sex, # 1 = male, 2 = female
    # change order of items, to make plot nicer later
    `Emotional Neglect`,
    `Physical Neglect`,
    `Emotional Abuse`,
    `Physical Abuse`,
    `Sexual Abuse`,
    everything()
  ) %>%
  bind_rows(biomarker_data_replication %>% # here we bind the replication sample to the original sample
              filter(MRID %in% relevant_IDs_replication$MRID) %>%
              select(.,
                     matches("RA1PRSEX|RA4QCT_EA|RA4QCT_SA|RA4QCT_PA|RA4QCT_EN|RA4QCT_PN|RA4Q4")) %>%
              `colnames<-`(c("Sex", var_names)) %>%
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
              select(Sex, # 1 = male, 2 = female
                     # change order of items, to make plot nicer later
                     `Emotional Neglect`,
                     `Physical Neglect`,
                     `Emotional Abuse`,
                     `Physical Abuse`,
                     `Sexual Abuse`,
                     everything()
              )) %>%  split(.$Sex) %>% 
 map(~select(., -c("Sex"))) # remove sex variable after grouping


# ----------------- estimate male network
graph_male <- estimateNetwork(
  graph_data_sex$`1`,
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

# ----------------- estimate female network
graph_female <- estimateNetwork(
  graph_data_sex$`2`,
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

# ----------------- run network comparison 
set.seed(1994)
compare_male_female <-
  NCT(
    graph_male,
    graph_female,
    it = 1000,
    test.edges = TRUE,
    edges = 'all',
    progressbar = TRUE
  )

p.adjust(compare_male_female$einv.pvals$`p-value`, method = "BH")
  
# ------------- Plotting networks
png(width=1250, height=450, "male_female_plot.png")
layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(2.5, 4))
qgraph(
  graph_male$graph,
  title.cex = 1.75,
  edge.width=1,
  layout = layout_network,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  legend=F,
  groups =  recode(
    wtc$membership,
    `2` = "Childhood Trauma",
    `1` = "Perceived Helplessness",
    `3` = "Perceived Self-Efficacy"
  ),
  legend.mode="style1",
  color =  c("grey",
             "#EBCC2A",
             "#78B7C5"
  ),
  label.cex = 1.45, 
  DoNotPlot=F,
  legend.cex = 0.69,
  vsize = 8,
  title = "Men",
  minimum = 0,
  maximum = 0.4814082
  
)

qgraph(
  graph_female$graph,
  layout = layout_network,
  edge.width=0.75,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  legend=T,
  GLratio = 2.25,
  groups =  recode(
    wtc$membership,
    `2` = "Childhood Trauma",
    `1` = "Perceived Helplessness",
    `3` = "Perceived Self-Efficacy"),
  layoutOffset = c(-0.305, 0),
  layoutScale = c(0.89, 1),
    legend.mode="style1",
    color =  c("grey",
               "#EBCC2A",
               "#78B7C5"
    ),
  label.cex = 1.375,
  vsize = 6.25,
  DoNotPlot=F,
    nodeNames = colnames(graph_original$graph),
  legend.cex = 0.65,
  title = "Women",
  minimum = 0,
  title.cex = 1.75,
  maximum = 0.4814082
  )

dev.off()


# ---------------------------- 5: Replication Analyses -----------------------------

# extract relevant variables from data set, basic "preprocessing" as above
graph_data_replication <- biomarker_data_replication %>%
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

# estimate replication network
graph_replication <- estimateNetwork(
  graph_data_replication,
  default = "ggmModSelect",
  tuning = 0,
  stepwise = TRUE,
  missing = "pairwise",
  corArgs = list(method = "spearman"),
  corMethod = "cor"
)

# basic comparison: correlate the two partial correlation matrices
cor(c(graph_original$graph), c(graph_replication$graph)) # 0.9242128

# Network Comparison for Differences in Structure, Global Strength & Individual Edges
set.seed(1337)
# may take longer to run on a standard PC!
compare_12 <-
  NCT(
    graph_original,
    graph_replication,
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

graph_replication_strength <-
  centralityTable(graph_replication$graph) %>% filter(measure == "Strength") %>%
  transmute(value)

cor(graph_original_strength, graph_replication_strength) # 0.9562717

# plot replication network & combined network
# estimate communities via walktrap
walktrap.community(as.igraph(qgraph(graph_replication$graph), attributes = TRUE))$membership

# estimate combined network
graph_combined <- estimateNetwork(
  rbind.data.frame(graph_data_original, graph_data_replication),
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
  layout = layout_network,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  legend=F,
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
             "#78B7C5"
  ),
  label.prop = 0.96,
  vsize=13,
  DoNotPlot=F,
  nodeNames = colnames(graph_original$graph),
  minimum = 0,
  maximum = 0.483,
  title="a) Original",
  title.cex=3
)

qgraph(
  graph_replication$graph,
  layout = layout_network,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  legend=F,
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
             "#78B7C5"
  ),
  label.prop = 0.96,
  vsize=13,
  DoNotPlot=F,
  nodeNames = colnames(graph_original$graph),
  minimum = 0,
  maximum = 0.483,
  title="b) Replication",
  title.cex=3
)

qgraph(
  graph_combined$graph,
  layout = layout_network,
  theme = "Borkulo",
  labels = c("EmN", "PhN", "EmA", "PhA", "SxA", 1:10),
  legend=F,
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
             "#78B7C5"
  ),
  label.prop = 0.96,
  vsize=14,
  DoNotPlot=F,
  nodeNames = colnames(graph_original$graph),
  minimum = 0,
  maximum = 0.483,
  title="c) Combined",
  title.cex=3
  
)

dev.off()
