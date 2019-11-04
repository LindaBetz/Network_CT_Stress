library(qgraph)
library(bootnet)
library(dplyr)


biomarker_data <- da29282.0001 


recode_vars <- c(
  "Confident to handle personal problems",
  # pos
  "Things were going your way",
  # pos
  "Able to control irritations in life",
  # pos
  "Felt on top of things"
  # pos
)
var_names <- c(
  "Upset by something unexpected",
  "Unable to control important things",
  "Felt nervous and stressed",
  "Confident to handle personal problems",
  # pos
  "Things were going your way",
  # pos
  "Could not cope with all things to do",
  "Able to control irritations in life",
  # pos
  "Felt on top of things",
  # pos
  "Angered by things outside control",
  "Difficulties piling up can't overcome",
  "Emotional Abuse",
  "Physical Abuse",
  "Sexual Abuse",
  "Emotional Neglect",
  "Physical Neglect"
)


graph_data <- biomarker_data %>%
  select(.,
         matches("B4QCT_EA|B4QCT_SA|B4QCT_PA|B4QCT_EN|B4QCT_PN|B4Q4")) %>%
  `colnames<-`(var_names) %>%
  mutate_all(as.numeric) 

 #mutate_at(recode_vars, ~recode(., `1` = 5, `2` = 4, `3` =3, `4` = 2, `5` = 1, .missing=NA_real_))

graph <- estimateNetwork(graph_data, 
                         default = "ggmModSelect", 
                         tuning = 0,
                         stepwise = TRUE,
                         missing = "pairwise",
                         corArgs=list(method="spearman"),
                         corMethod = "cor")

graph_plot <- qgraph(graph$graph,
       layout = lay$layout,
       theme="Borkulo",
       labels = TRUE,
       labels = c(1:15),
       cut=0,
       GLratio = 1.9,
       layoutOffset = c(-0.26, 0),
       layoutScale = c(0.9, 1),
       nodeNames = var_names,
       legend.cex=0.49,
       legend = TRUE)

set.seed(1)
g <- igraph::as.igraph(graph_plot, attributes=TRUE)
groups <- igraph::spinglass.community(g)$membership %>% as.character()

graph_plot <- qgraph(graph$graph,
                     layout = "spring",
                     theme="Borkulo",
                     labels = TRUE,
                     layout.par = list(repulse.rad = 500, max.delta = 15),
                     labels = c(1:15),
                     GLratio = 1.9,
                     groups=groups,
                     legend.mode = "style2",
                     layoutOffset = c(-0.26, 0),
                     layoutScale = c(0.9, 1),
                     palette = "colorblind",
                     nodeNames = var_names,
                     legend.cex=0.49,
                     legend = TRUE)

results_case <- bootnet(graph, type="case", nBoots=500)
results_boot <- bootnet(graph, nBoots=500)
