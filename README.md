# Transcriptomic-data-analysis-DEMO
## Overview
It is a common strategy to analyze transcriptomic data from samples subjected to different treatments (e.g., resistant vs. susceptible) in order to identify gene expression profiles, discover key genes, and infer related biological pathways. Here, a straightforward case of mining acaricide resistance-related genes and pathways using the transcriptomic data of _Tetranychus urticae_ is presented.

## Workflow
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/flowchart.png" width="450">

Data in **bold** format are acquired from online database or made by myself following the format.

Data source:  
Transcriptomic data: From NCBI (Accession: PRJNA946758).  
Reference genome: From https://figshare.com/articles/dataset/Chromosome-level_genome_assembly_of_the_two-spotted_spider_mite_i_Tetranychus_urticae_i_/25241794/1.  
GFF file: Same as Reference genome.  
Desination table: Self-made .csv or .txt file, following the format of "Experiment design file" in iDEP.  

In this demo, I use three packages/web tools to improve the workflow of the analysis:  
1. Trimming sequence and generate QC report by _**AfterQC**_ due to its versatility in quality control, data filtering, error profiling and base correction automatically.  
2. Annotating sequence by _**eggNOG-mapper**_ which is a powerful function annotation and prediction web tool.  
3. Analyzing the gene counts matrix by _**iDEP**_, a web application with easy-to-use GUI for differential expression and pathway analysis, which saving a lot of hours to understand R codes and debug.  

A tutorial with example files is available in Tutorial folder.  

## Results
(1) Preprocessing:  
Little bias was observed in the raw read counts, and the distribution of rlog-transformed data was consistent across samples and treatments.    

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/raw_counts_barplot.png" width="450">

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/transformed_boxplot.png" width="450">

(2) Clustering analysis:  
The expression patterns of "WITHEXPOSURE" (WE) samples were similar to those of "WITHOUTEXPOSURE" (WOE), while "REFERENCE" (REF) exhibited a distinct difference compared to the others.  

Hierarchical clustering    |  K-means clustering
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/heatmap_main_H.png"  height=250>|<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/heatmap_main_K.png" height=250>

(3) Principal component analysis (PCA)  
In the results of the PCA, WE and WOE were clustered on the right side, while REF was positioned on the left side. This arrangement revealed a distinct pattern, indicating that REF was different from WE and WOE, whereas a subtle difference existed between WE and WOE.

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/pca_plot.png" width="450">

(4) Differential expression analysis  
In the comparison between WE and REF, 781 up-regulated genes and 500 down-regulated genes were identified based on statistical thresholds (log fold change ≥ 2 and FDR ≤ 0.05). In contrast, the comparison of WOE and REF revealed 787 up-regulated genes and 500 down-regulated genes using the same criteria.   

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/sig_gene_stats.png" width="450">

The volcano plots for the two comparisons were shown below. The distribution of genes in both comparisons was similar, indicating that there is little difference between WE and WOE.

WITHEXPOSURE vs REFERENCE  |  WiITHOUTEXPOSURE vs REFERENCE
:-------------------------:|:-------------------------:
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/volcano_plot_WE_REF.png" width="250">|<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/volcano_plot_WOE_REF.png" width="250">

According to the analysis of the top five up-regulated and down-regulated genes in both comparisons, the functions of the up-regulated genes primarily involve DUF4774 domain-containing proteins, enoyl-reductases, protein ubiquitination, and C2H2-type domain-containing proteins. In contrast, the functions of the down-regulated genes are mainly associated with the lipocalin/cytosolic fatty-acid binding protein family, transferase activity, and the negative regulation of lipoprotein oxidation.

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/top5function_for_updown.png" width="450">

(5) Over-representation analysis  
The up-regulated pathways are primarily observed in glutathione metabolic processes across both comparisons, indicating an enhanced cellular response to oxidative stress. Notably, the "WE vs. REF" comparison shows up-regulation of UDP-glycosyltransferase activity, monooxygenase activity, and lysosome pathways. These findings suggest a heightened defense response and detoxification processes triggered by exposure. In the "WOE vs. REF" group, pathways related to proteolysis, cellular modified amino acid metabolism process, and sphingolipid metabolism process are also up-regulated, highlighting a potential shift in cellular energy and lipid signaling dynamics in the absence of external exposure.

Conversely, the down-regulated pathways encompass essential cellular and molecular components, including integral and intrinsic component of membrane, iron ion binding, and catalytic activity, particularly in the "WE vs. REF" group. Processes such as the TRAMP complex and peptidoglycan muralytic activity are significantly reduced in the samples of WE, indicating a suppressive effect on RNA and peptidoglycan.  

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/pathway_ora.png" width="450">

Note: NPA is the abbrevation for No Pathway Available.

(6) Pathway analysis  
 In the up-regulated pathways, nucleosome organization, phosphatase regulation, and glycogen biosynthetic processes are significantly activated in the "WE vs. REF", suggesting an increase in chromatin modification and metabolic activity. Molecular functions such as oxidoreductase activity and aspartic-type endopeptidase activity are also up-regulated, reflecting heightened enzymatic activity. Additionally, KEGG pathways related to phenylalanine metabolism and drug metabolism by cytochrome P450  are prominent. The "WOE vs. REF" comparison underscores the activation of positive TOR signaling regulation and arachidonic acid secretion.  

In terms of down-regulated pathways, "WE vs. REF" reveals reductions in processes such as cell junction assembly, RNA polyadenylation, and miRNA-mediated gene silencing. This indicates a decrease in cellular structural organization and gene silencing mechanisms. Molecular functions like intramolecular transferase and mannosyltransferase activity, along with KEGG pathways related to N-glycan biosynthesis and mismatch repair, are also down-regulated. Similarly, in the "WOE vs. REF" group, the down-regulation occurs in pathways such as cell junction organization, RNA polyadenylation, and various types of glycan biosynthesis, indicating less structural and metabolic activity in this condition.  

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Content/pathway_gsea.png" width="450">

Note: NPA is the abbrevation for No Pathway Available.

## Reference
Anders, S., Pyl, P. T., & Huber, W. (2015). HTSeq—a Python framework to work with high-throughput sequencing data. _Bioinformatics_, 31(2), 166-169.  
Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. _Molecular Biology and Evolution_, 38(12), 5825-5829.  
Chen, S., Huang, T., Zhou, Y., Han, Y., Xu, M., & Gu, J. (2017). AfterQC: automatic filtering, trimming, error removing and quality control for fastq data. _BMC Bioinformatics_, 18, 91-100.  
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., & Davies, R. M. (2021). Twelve years of SAMtools and BCFtools. _Gigascience_, 10(2), giab008.  
Ge, S. X., Son, E. W., & Yao, R. (2018). iDEP: an integrated web application for differential expression and pathway analysis of RNA-Seq data. _BMC Bioinformatics_, 19, 1-24.  
Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. _Nature Methods_, 12(4), 357-360.  
Pertea, G., & Pertea, M. (2020). GFF utilities: GffRead and GffCompare. _F1000Research_, 9.  
SRA Toolkit Development Team. SRA Toolkit. https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software.  
Wei, Shu-Jun; Cao, Li-Jun (2024). Chromosome-level genome assembly of the two-spotted spider mite _Tetranychus urticae_. _figshare_. Dataset. https://doi.org/10.6084/m9.figshare.25241794.v1.  
