library(qgraph)
library(bootnet)
library(dplyr)


biomarker_data <- da29282.0001 



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
  #mutate_at(recode_vars, ~recode(., `1` = 5, `2` = 4, `3` =3, `4` = 2, `5` = 1))


graph <- estimateNetwork(graph_data, 
                         default = "ggmModSelect", 
                         tuning = 0,
                         corMethod = "cor")

qgraph(graph$graph,
       layout = "spring",
       theme="Borkulo",
       labels = TRUE,
       labels = c(1:15),
       threshold = 0.05,
       cut=0,
       GLratio = 1.9,
       legend.mode = "style2",
       layoutOffset = c(-0.1, 0),
       layoutScale = c(0.9, 0.9),
       palette = "colorblind",
       groups = c(rep("1) Perceived Stress", 10), rep("2) Childhood Trauma", 5)),
       nodeNames = var_names,
       legend.cex=0.55,
       legend = TRUE)
