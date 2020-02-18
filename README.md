# Network of Childhood Trauma & Perceived Stress

R-code to reproduce analyses described in 
"Relationships between childhood trauma and subjective experiences of stress in the general population: a network perspective"
by Linda T. Betz, Nora Penzel, Marlene Rosen, Joseph Kambeitz

Code by L. Betz (linda.betz@uk-koeln.de)

There are three files: 

1) Code_Main_Analysis.R provides code to reproduce findings and plots reported in the main manuscript.

2) Code_Supplementary_Analysis.R provides code to reproduce findings and plots reported in the supplementary file.

3) Weighted_adjacency_matrices.zip provides the weighted adjacency matrices containing the partial correlation coefficients for the Gaussian Graphical Models reported in the paper as .csv-files.


Data required for the analysis is available at https://www.icpsr.umich.edu/icpsrweb/:

MIDUS Biomarker: https://doi.org/10.3886/ICPSR29282.v9

MIDUS Biomarker Refresher: https://doi.org/10.3886/ICPSR36901.v6

MIDUS II: https://doi.org/10.3886/ICPSR04652.v7

MIDUS Refresher: https://doi.org/10.3886/ICPSR36532.v3


For optimal reproducibility, we use the R-package "checkpoint", snapshot date November 5, 2019 (R version 3.6.1). Information on how this works can be found at https://mran.microsoft.com/documents/rro/reproducibility.
