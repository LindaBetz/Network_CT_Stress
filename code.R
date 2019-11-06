

# ---------------------------- 1: Load libraries & packages -----------------------------
library(ggsci)
library(qgraph)
library(bootnet)
library(purrr)
library(dplyr)

# all data sets are available at https://www.icpsr.umich.edu/icpsrweb/

biomarker_data <- da29282.0001
wave_2_data <- da04652.0001
biomarker_refresher_data <- da36901.0001
wave_2_refresher_data <- 

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
relevant_IDs <- biomarker_data %>%
  select(.,
         matches("M2ID|B4QCT_EA|B4QCT_SA|B4QCT_PA|B4QCT_EN|B4QCT_PN|B4Q4")) %>%
  mutate(na_per_row = rowSums(is.na(.) / 15)) %>%
  filter(na_per_row != 1) %>%
  transmute(M2ID) # 3 people have missing data in all variables of interest

# create data set to calculate sample statstiscs
desc_data <- biomarker_data %>%
  left_join(wave_2_data, by = "M2ID") %>%
  filter(M2ID %in% relevant_IDs$M2ID) %>%
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
  )


# frequency/count data
desc_data %>%
  select(., c(Site, MIDUS_Sample, Ethnicity, Education)) %>%
  map( ~ table(.) / sum(!is.na(.))) %>%
  map( ~ round(., 3))

# continous data
desc_data %>%
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
graph_data_original <- biomarker_data %>%
  # filter(SAMPLMAJ == "(01) MAIN RDD") %>%
  select(.,
         matches("B4QCT_EA|B4QCT_SA|B4QCT_PA|B4QCT_EN|B4QCT_PN|B4Q4")) %>%
  `colnames<-`(var_names) %>%
  mutate_all(as.numeric) %>%
  mutate_at(recode_vars,
            ~ recode( # recode positive items
              .,
              `1` = 5,
              `2` = 4,
              `3` = 3,
              `4` = 2,
              `5` = 1,
              .missing = NA_real_
            )) %>%
  select( # change order of items, to make plot nicer later
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

# plot network
graph_plot <- qgraph(
  graph$graph,
  layout = "spring",
  theme = "Borkulo",
  labels = TRUE,
  layout.par = list(max.delta = 10),
  labels = c(1:15),
  GLratio = 1.9,
  groups = c(rep("1) Childhood Trauma", 5), rep("2) Perceived Stress", 10)),
  layoutOffset = c(-0.25, 0),
  layoutScale = c(0.9, 1),
  label.cex = 1.25,
  filetype = "png",
  color = c(rep("#D3D3D3", 5),
            rep("#5DBCD2", 10)),
  nodeNames = colnames(graph_data),
  legend.cex = 0.67,
  legend = TRUE
)

# ---------------------------- 4: Sensitvity Analyses -----------------------------
# case-drop bootstrapping
results_case <- bootnet(graph, type = "case", nBoots = 1000)
corStability(results_case) # 0.75 for edge & strength
plot(results_case)

# regular bootstrapping for edge weights
results_boot <- bootnet(graph, nBoots = 1000)
plot(results_boot, order = "sample", split0 = TRUE)


# ---------------------------- 5: Replication Analyses -----------------------------
# extract relevant variables from data set, basic "preprocessing" as above
graph_data_refresh <- biomarker_refresher_data %>%
  # filter(SAMPLMAJ == "(01) MAIN RDD") %>%
  select(.,
         matches("RA4QCT_EA|RA4QCT_SA|RA4QCT_PA|RA4QCT_EN|RA4QCT_PN|RA4Q4")) %>%
  `colnames<-`(var_names) %>%
  mutate_all(as.numeric) %>%
  mutate_at(recode_vars,
            ~ recode(  # recode positive items
              .,
              `1` = 5,
              `2` = 4,
              `3` = 3,
              `4` = 2,
              `5` = 1,
              .missing = NA_real_
            )) %>%
  select( # change order of items, to make plot nicer later
    `Emotional Neglect`,
    `Physical Neglect`,
    `Emotional Abuse`,
    `Physical Abuse`,
    `Sexual Abuse`,
    everything()
  )

# estimate network
graph_refresher <- estimateNetwork(
  graph_data,
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