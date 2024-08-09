# RNA-seq-data-analysis-DEMO
## Overview
It's a common strategy to analyze RNA-seq data of samples with different treatments (_e.g._ Resistance vs Susceptible) to identify gene expression profile, find key genes and biological pathways. Here, a simple case of mining resistance-related genes and pathways using RNA-seq data of _Tetranychus urticae_ is presented. 

## Workflow
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/flowchart.png" width="400">

## Results
(1). "With exposure" samples is similar to "Without exposure" ones, while "Reference" ones show different pattern compared to the others (see Figiure 1, 2 and 3).  
(2). Up-regulated differential expressed genes are mainly associated to molecular transportation, detoxification, ubiquitination and hydrolase/proteolysis, while down-regulated ones are associated to ATPase, peptidase inhibitor and reduction of peroxides (see Table 1).  
(3). Up-regulated pathways are mainly associated to membrane of cells/organelles (cellular component), glucan/ATP biosynthetic process (biological process), UDP-glycosyltransferase/monooxygenase/aspartic-type-peptidase activity (molecular function) and drug/glutathione/pyruavte metabolism (KEGG); down-regulated pathways  are mainly associated to cytoplasma/Golgi stack  (cellular component), mannosylation/DNA-topological-change (biological process), acyl-CoA/pigment binding (molecular function) and N-glycan-metabolism/autophagy (KEGG) (see Table 2).

### Figure 1. Heatmap
Hierarchical clustering    |  K-means clustering
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/heatmap_H.png"  height=250>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/heatmap_K.png" height=250>

### Figure 2. Differential expressed genes bar chart
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/sig_gene_stats.png"  height=250>

### Figure 3. Volcano plot
With Exposure vs Reference  |  Without Exposure vs Reference
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_we_ref.png" height=250>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_woe_ref.png" height=250>

### Table 1. Top 5 up/down regulated differential expressed genes for each comparison
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/top5function_for_updown.png" height=250>
Ranked by log(Fold change).

### Table 2. Up/down regulated pathways for each comparison
Over-representation analysis    |  Gene set enrichment analysis
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/pathway_ora.png" height=150>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/pathway_gesa.png" height=150>


## Details of each step of flow chart
The design of RNA-seq data contains 3 treatments, including "With exposure", "Without exposure" and "Reference". The procedure of analysis consist of 5 steps: (1) Download from database / custom data, (2) Quality control, (3) Reference genome indexing, (4) Mapping to genome, (5) Read counts calculation and (6) DEG & enrichment analysis.

## Note
In this demo, I use three packages/web tools to improve the workflow of the analysis: 
1. Trimming sequence and generate QC report by **AfterQC** due to its versatility in quality control,
 data filtering, error profiling and base correction automatically. 
2. Annotating sequence by **eggNOG-mapper** which is a powerful function annotation and prediction web tool. 
3. Analyzing the gene counts matrix by **iDEP**, a web application with easy-to-use GUI for differential expression and pathway analysis, which saving a lot of time to understand R codes and debug.

## Reference
