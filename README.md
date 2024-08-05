# RNA-seq-data-analysis-DEMO
## Overview
It's a common strategy to analyze RNA-seq data of samples with different treatments (_e.g._ Resistance vs Susceptible) to identify gene expression profile, find key genes and biological pathways. Here, a simple case of mining resistance-related genes using RNA-seq data of _Tetranychus urticae_ is presented. 

## Workflow
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/flowchart.png" width="400">

## Results



## Details of each step of flow chart
The design of RNA-seq data contains 3 treatments, including "With exposure", "Without exposure" and "Reference". The procedure of analysis consist of 5 steps: (1) Download from database / custom data, (2) Quality control, (3) Reference genome indexing, (4) Mapping to genome, (5) Read counts calculation and (6) DEG & enrichment analysis.

## Difference between my strategy and common ones
1. Trimming sequence and generate QC report by **AfterQC** due to its versatility in quality control,
 data filtering, error profiling and base correction automatically. 
2. Annotating sequence by **eggNOG-mapper** which is a powerful function annotation and prediction web tool. 
3. Analyzing the gene counts matrix by **iDEP**, a web application with easy-to-use GUI for differential expression and pathway analysis, which saving a lot of time to understand R codes and debug.

## Figures and Tables
### Heatmap
Hierarchical clustering    |  K-means clustering
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/heatmap_H.png"  height=250>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/heatmap_K.png" height=250>

### Differential expressed genes bar chart
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/sig_gene_stats.png"  width=300>

### Volcano plot
Without Exposure vs Reference   |  With Exposure vs Without Exposure   |  With Exposure vs Reference 
:-------------------------:|:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_woe_ref.png">|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_we_woe.png">|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_we_ref.png">

### Top 5 up/down regulated differential expressed genes for each comparison
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/top5function_for_updown.png" width=300>
Filtered by log(Fold change)

### Up/down regulated pathways for each comparison
Over-representation analysis    |  Gene set enrichment analysis
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/pathway_ora.png" height=150>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/pathway_gesa.png" height=150>
