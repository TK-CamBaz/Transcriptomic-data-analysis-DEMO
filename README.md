# RNA-seq-data-analysis-DEMO
## Overview
A hands-on work of RNA-seq data analysis for mining resistance-related genes in agriculture pests is presented. The purpose is to explore potential genes related to acaricide resistance in "Tetranychus urticae" using RNA-seq data and bioinformatic methods. The design of RNA-seq data contains 3 treatments, including "With exposure", "Without exposure" and "Reference". The procedure of analysis consist of 5 steps: (1) Download from database / custom data, (2) Quality control, (3) Reference genome indexing, (4) Mapping to genome, (5) Read counts calculation and (6) DEG & enrichment analysis.

## Difference between my strategy and common ones
1. Trimming sequence and generate QC report by **AfterQC** due to its versatility in quality control,
 data filtering, error profiling and base correction automatically. 
2. Annotating sequence by **eggNOG-mapper** which is a powerful function annotation and prediction web tool. 
3. Analyzing the gene counts matrix by **iDEP**, a web application with easy-to-use GUI for differential expression and pathway analysis, which saving a lot of time to understand R codes and debug.

## Flow chart
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/flowchart.png" width="400">

## Results
### Figures and Tables
Hierarchical clustering    |  K-means clustering
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/heatmap_H.png"  height=250>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/heatmap_K.png" height=250>

<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/sig_gene_stats.png"  width=300>


