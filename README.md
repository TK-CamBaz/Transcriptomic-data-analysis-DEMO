# Transcriptomic-data-analysis-DEMO
## Overview
It is a common strategy to analyze transcriptomic data from samples subjected to different treatments (e.g., resistant vs. susceptible) in order to identify gene expression profiles, discover key genes, and infer related biological pathways. Here, a straightforward case of mining acaricide resistance-related genes and pathways using the transcriptomic data of Tetranychus urticae is presented.

## Workflow
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/flowchart.png" width="450">

Data in **bold** format are acquired from online database or made following the format.

Data source:  
Transcriptomic data: From NCBI (Accession: PRJNA946758).  
Reference genome: From https://figshare.com/articles/dataset/Chromosome-level_genome_assembly_of_the_two-spotted_spider_mite_i_Tetranychus_urticae_i_/25241794/1.  
GFF file: Same as Reference genome.  
Desination table: Self-made .csv or .txt file, following the format of "Experiment design file" in iDEP.

## Results
(1) Preprocessing:  
Little bias was observed in the raw read counts, and the distribution of rlog-transformed data was consistent across samples and treatments.    

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/raw_counts_barplot.png" width="450">

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/transformed_boxplot.png" width="450">

(2) Clustering analysis:  
The expression patterns of "WITHEXPOSURE" (WE) samples were similar to those of "WITHOUTEXPOSURE" (WOE), while "REFERENCE" (REF) exhibited a distinct difference compared to the others.  

Hierarchical clustering    |  K-means clustering
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/heatmap_main_H.png"  height=250>|<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/heatmap_main_K.png" height=250>

(3) Principal component analysis (PCA)  
In the results of the PCA, WE and WOE were clustered on the right side, while REF was positioned on the left side. This arrangement revealed a distinct pattern, indicating that REF was different from WE and WOE, whereas a subtle difference existed between WE and WOE.

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/pca_plot.png" width="450">

(4) Differential expression analysis  
In the comparison between WE and REF, 781 up-regulated genes and 500 down-regulated genes were identified based on statistical thresholds (log fold change ≥ 2 and FDR ≤ 0.05). In contrast, the comparison of WOE and REF revealed 787 up-regulated genes and 500 down-regulated genes using the same criteria.   

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/sig_gene_stats.png" width="450">

The volcano plots of the two comparison were shown below. 

WITHEXPOSURE vs REFERENCE  |  WiITHOUTEXPOSURE vs REFERENCE
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/volcano_plot_WE_REF.png" width="250">|<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Contents/volcano_plot_WOE_REF.png" width="250">

(2). Little difference is observed between WE_Ref and WOE_Ref, wihch also suggests the similarity between WE and WOE (see Figiure 2, 3 and Table 1).  
(3). Up-regulated differential expressed genes are mainly associated to molecular transportation, detoxification, ubiquitination and hydrolase/proteolysis, while down-regulated ones are associated to ATPase, peptidase inhibitor and reduction of peroxides (see Table 1).  
(4). Up-regulated pathways are mainly associated to membrane of cells/organelles (cellular component), glucan/ATP biosynthetic process (biological process), UDP-glycosyltransferase/monooxygenase/aspartic-type-peptidase activity (molecular function) and drug/glutathione/cytochrome-P450 metabolism (KEGG); down-regulated pathways are mainly associated to cytoplasma/Golgi stack  (cellular component), DNA topological change (biological process), various-transferase/DNA-topoisomerase (molecular function) and N-glycan-metabolism/autophagy (KEGG) (see Table 2).  

### Figure 1. Heatmap


### Figure 2. Differential expressed genes bar chart
Conditions: FDR <= 0.05 & Log(FC) >=2.  

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/FigureTable/sig_gene_stats.png"  height=250>

### Figure 3. Volcano plot
Conditions: FDR <= 0.05 & Log(FC) >=2.  

With Exposure vs Reference  |  Without Exposure vs Reference
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_we_ref.png" height=200>|<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_woe_ref.png" height=200>

### Table 1. Top 5 up/down regulated differential expressed genes for each comparison
Conditions: Ranked by log(FC).  

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/FigureTable/top5function_for_updown.png" height=250>

### Table 2. Up/down regulated pathways for each comparison
Conditions: Ranked by FDR.

Over-representation analysis    |  Gene set enrichment analysis
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/FigureTable/pathway_ora.png" height=150>|<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/FigureTable/pathway_gesa.png" height=150>

## Note
In this demo, I use three packages/web tools to improve the workflow of the analysis: 
1. Trimming sequence and generate QC report by **AfterQC** due to its versatility in quality control,
 data filtering, error profiling and base correction automatically. 
2. Annotating sequence by **eggNOG-mapper** which is a powerful function annotation and prediction web tool. 
3. Analyzing the gene counts matrix by **iDEP**, a web application with easy-to-use GUI for differential expression and pathway analysis, which saving a lot of time to understand R codes and debug.

## Reference
Anders, S., Pyl, P. T., & Huber, W. (2015). HTSeq—a Python framework to work with high-throughput sequencing data. Bioinformatics, 31(2), 166-169.  
Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. _Molecular Biology and Evolution_, 38(12), 5825-5829.  
Chen, S., Huang, T., Zhou, Y., Han, Y., Xu, M., & Gu, J. (2017). AfterQC: automatic filtering, trimming, error removing and quality control for fastq data. _BMC Bioinformatics_, 18, 91-100.  
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., & Davies, R. M. (2021). Twelve years of SAMtools and BCFtools. _Gigascience_, 10(2), giab008.  
Ge, S. X., Son, E. W., & Yao, R. (2018). iDEP: an integrated web application for differential expression and pathway analysis of RNA-Seq data. _BMC Bioinformatics_, 19, 1-24.  
Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. _Nature Methods_, 12(4), 357-360.  
Pertea, G., & Pertea, M. (2020). GFF utilities: GffRead and GffCompare. _F1000Research_, 9.  
SRA Toolkit Development Team. SRA Toolkit. https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software.  
Wei, Shu-Jun; Cao, Li-Jun (2024). Chromosome-level genome assembly of the two-spotted spider mite _Tetranychus urticae_. _figshare_. Dataset. https://doi.org/10.6084/m9.figshare.25241794.v1.  
