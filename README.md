# RNA-seq-data-analysis-DEMO
## Overview
It's a common strategy to analyze RNA-seq data of samples with different treatments (_e.g._ Resistance vs Susceptible) to identify gene expression profile, find key genes and biological pathways. Here, a simple case of mining resistance-related genes and pathways using RNA-seq data of _Tetranychus urticae_ is presented. 

## Workflow
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/flowchart.png" width="450">

Data source:  
RNAseq data: From NCBI (Accession: PRJNA946758).  
Reference genome: From https://figshare.com/articles/dataset/Chromosome-level_genome_assembly_of_the_two-spotted_spider_mite_i_Tetranychus_urticae_i_/25241794/1.  
GFF file: Same as Reference genome.  
Desination table: Self-made .csv or .txt file, following the format of "Experiment design file" in iDEP.

## Results
(1). Expression patterns of "With exposure" (WE) samples is similar to ones of "Without exposure" (WOE) , while "Reference" (Ref) shows distinct difference compared to the others (see Figiure 1).  
(2). Little difference is observed between WE_Ref and WOE_Ref, wihch also suggests the similarity between WE and WOE (see Figiure 2, 3 and Table 1).  
(3). Up-regulated differential expressed genes are mainly associated to molecular transportation, detoxification, ubiquitination and hydrolase/proteolysis, while down-regulated ones are associated to ATPase, peptidase inhibitor and reduction of peroxides (see Table 1).  
(4). Up-regulated pathways are mainly associated to membrane of cells/organelles (cellular component), glucan/ATP biosynthetic process (biological process), UDP-glycosyltransferase/monooxygenase/aspartic-type-peptidase activity (molecular function) and drug/glutathione/cytochrome-P450 metabolism (KEGG); down-regulated pathways are mainly associated to cytoplasma/Golgi stack  (cellular component), DNA topological change (biological process), various-transferase/DNA-topoisomerase (molecular function) and N-glycan-metabolism/autophagy (KEGG) (see Table 2).  

### Figure 1. Heatmap
Hierarchical clustering    |  K-means clustering
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/heatmap_H.png"  height=250>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/heatmap_K.png" height=250>

### Figure 2. Differential expressed genes bar chart
Conditions: FDR <= 0.05 & Log(FC) >=2.  

<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/sig_gene_stats.png"  height=250>

### Figure 3. Volcano plot
Conditions: FDR <= 0.05 & Log(FC) >=2.  

With Exposure vs Reference  |  Without Exposure vs Reference
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_we_ref.png" height=200>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/volcano_plot_woe_ref.png" height=200>

### Table 1. Top 5 up/down regulated differential expressed genes for each comparison
Conditions: Ranked by log(FC).  

<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/top5function_for_updown.png" height=250>

### Table 2. Up/down regulated pathways for each comparison
Conditions: Ranked by FDR.

Over-representation analysis    |  Gene set enrichment analysis
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/pathway_ora.png" height=150>|<img src="https://github.com/TK-CamBaz/RNA-seq-data-analysis-DEMO/blob/main/FigureTable/pathway_gesa.png" height=150>

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
