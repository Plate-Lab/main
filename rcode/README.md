# R Code
This folder contains the Plate Lab's commonly-used R code. We currently have two pipelines in R for the analysis and visualization of both DDA/TMT and DIA LC-MS/MS data.

## Pipelines
 - [DIA](DIA.R)
 - [DDA TMT](DDA_TMT.R)

### General Workflow
1. Load data file
2. Run QC - q.value, contaminant removal
3. Normalize - median normalize and/or bait normalize, log2 transformation
4. Statistics - ANOVA and/or Volcano Plots
5. Calculate fold changes and/or avg.FC per protein
6. K-means clustering, heatmap
