# Analysis Workflow 
## 1. Data Preparation
In this project, I have prepared three files: transcriptome data, whole-genome sequence, and genome annotation. The transcriptome data is from NCBI (accession: PRJNA946758), and the whole-genome sequence and annotation are from Figshare (https://figshare.com/articles/dataset/Chromosome-level_genome_assembly_of_the_two-spotted_spider_mite_i_Tetranychus_urticae_i_/25241794/1). The transcriptome data contains 12 SRA files, each containing transcriptome data for different samples, sequenced using the Illumina NovaSeq 6000 platform with paired-end sequencing layout. The experimental design includes one control group (susceptible Tetranychus urticae strain not exposed to acaricides) and two treatment groups (resistant Tetranychus urticae strain continuously exposed to acaricides, resistant strain not exposed to acaricides), each with four replicates. The susceptible strain is the London reference strain, and the resistant strain is FP9 (highly resistant to acequinocyl and bifenazate). The acaricide exposure is acequinocyl. The whole-genome sequence is the chromosome-level assembly of Tetranychus urticae, provided as a FASTA file. The genome annotation provides the coding sequences (CDS) of Tetranychus urticae on the genome, provided as a GFF file.  
 The sra files can be downloaded using the prefetch command from the sra-toolkit:  
For example, with SRX19731018: 
```
prefetch SRX19731018
```
The whole genome sequence and annotations can be downloaded directly from the provided link.

## 2. Data Preprocessing
### (1) Extract FASTQ files from SRA files and perform quality control
First, use parallel-fastq-dump from the SRA toolkit to extract FASTQ files from SRA files:  
For example, with SRX19731018.sra:
```
parallel-fastq-dump -s SRX19731018.sra -t 8 --split-3
```
In the case of paired-end sequencing, two FASTQ files will be generated, both of which should be retained.

Next, use AfterQC to clean the sequences, mainly to remove adapters and sequences that are too short. AfterQC automatically determines the cleaning conditions and outputs reports comparing the sequences before and after cleaning.  
For example, with SRR23919699:
```
cd ~/AfterQC-0.9.7
pypy after.py -1 SRR23919699_1.fastq -2 SRR23919699_2.fastq
```
The output includes three directories: good, bad, and QC. The fq files in the good folder are the clean reads, the fq files in the bad folder are the removed reads, and the HTML files in the QC folder are the reports that compare the sequences before and after cleaning.

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/AfterQC_partial_res.png" width="450">

### (2) Build a whole genome index
The purpose of this step is to speed up the sequence mapping process, which is done using HISAT2.  
For example, with tu.genome.ipm_v2.fasta:
```
hisat2-build -p 8 tu.genome.ipm_v2.fasta index/tu.genome_index
```
The output is a directory named tu.genome_index, containing eight H2 files.

### (3) Mapping sequences to the reference genome
Use HISAT2 for this step.  
For example, with SRR23919699_1.good.fq and SRR23919699_2.good.fq:
```
hisat2 -t -p 2 -x index/tu.genome_index -1 SRR23919699_1.good.fq -2 SRR23919699_2.good.fq -S SRR23919699.sam
```
The output is a SAM file.

Next, convert the SAM file to a BAM file. This step is a format conversion that significantly reduces the file size without altering its content.
Use samtools for this step.  
For example, with SRR23919699.sam:
```
samtools sort -@2 -m 200M -o SRR23919699.bam SRR23919699.sam
```
The output is a BAM file.

Then, use qualimap to check the mapping results:
```
cd ~/qualimap_v2.3
./qualimap
```
<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/qualimap_GUI.png" width="450">

After the window opens, select the BAM file for analysis. Once complete, a report will be displayed, which can also be exported as a PDF or HTML.

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/qualimap_res.png" width="450">

You can also generate a BAI file from the BAM file, which can be loaded into IGV (Integrative Genomics Viewer) along with the BAM file for visualizing the mapping results: 
```
samtools index SRR23919699.bam SRR23919699.bai
```
Then open IGV:
```
igv
```

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/igv_GUI.png" width="450">

Select the BAM file and the whole genome sequence. Make sure the BAM file and BAI file are in the same directory.

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/igv_res.png" width="450">

### (4) Quantify the transcripts in each sample after sequence mapping
This step requires the BAM file and a GTF file. First, convert the whole genome annotation GFF file to a gtf file using gffread:
For example, with tu_evm_out.gff3:
```
gffread tu_evm_out.gff3 -T -o tu_evm_out.gtf
```
Then, use HTSeq for counting the reads:  
For example, with SRR23919699.bam:  
```
htseq-count -q -f bam -s no -i gene_id SRR23919699.bam tu_evm_out.gtf > SRR23919699.count
```
Repeat the above steps for all samples to generate count files. Use Excel or other software (R or Python, etc.) to merge the files.

### (5) Obtain the count matrix
After merging the count files, the first column in the resulting file contains transcript IDs, which are not valid Gene IDs. Therefore, use an online database to convert the transcript IDs to Gene IDs. First, extract the sequences for each transcript using the whole genome sequence and whole genome annotation (in gtf format):
```
gffread -w transcripts.fa -g tu.genome.ipm_v2.fasta tu_evm_out.gtf
```
Then, use eggNOG-mapper (http://eggnog-mapper.embl.de/), select CDS, upload the transcript sequence file (transcripts.fa), provide your email, and click Submit. Once the comparison is complete, download the csv file (select csv; although the downloaded file has a tsv extension, the content is not affected). Then, use Merge.annotation.counts.R (or Excel's vlookup function) to generate the count matrix file.

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/emapper_UI.png" width="450">

## 3. Analyze the count matrix using iDEP
iDEP is a web tool designed for RNA-seq data analysis, offering common analysis methods such as data exploration, differential gene expression analysis, over-representation analysis, and gene set enrichment analysis. Besides uploading the count matrix, you can also prepare and upload an experimental design file. The format can be checked by clicking Info. In this experiment, the material used is the two-spotted spider mite, so make sure to select the correct species in the Species section. 

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_UI.png" width="450">

### (1) Preprocessing
Click Pre-Process, and each button allows you to view the distribution of Read counts across samples. Selecting rlog for expression value transformation helps to improve the correlation among samples under the same treatment. 

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_preprocess.png" width="450">

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_preprocess_2.png" width="450">

### (2) Clustering analysis
Click Clustering. Two options are available: hierarchical clustering and k-means clustering. The clustering results are presented as a heatmap showing the gene expression profiles of each sample. You can select regions of interest on the heatmap to generate a sub-heatmap, which allows further investigation of the expression profiles. The sub-heatmap also lists gene names for easy tracking and recording.

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_clustering.png" width="450">

### (3) Principal component analysis
Click PCA to observe the similarity between samples. Other dimensionality reduction methods are also available.  

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_Dimension_reduction.png" width="450">

### (4) Differential gene expression analysis
Click DEG1, select the comparison groups, and click Submit. The results are displayed as bar plots and tables showing the number of upregulated and downregulated differentially expressed genes in each comparison. The results can be downloaded from the Result & data section. 

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_DEG.png" width="450">

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_DEG_2.png" width="450">

### (5) Differential gene expression visualization and over-representation analysis
Click DEG2. You can observe the distribution of differentially expressed genes across the genome using heatmaps, volcano plots, MA plots, and scatter plots. 

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_DEG_3.png" width="450">

Select Enrichment to perform over-representation analysis. On the left, choose the group to analyze. The differentially expressed genes in this group will be treated as a gene list and compared to the selected gene sets from the database using a hypergeometric test. Pathways with an adjusted p-value (FDR) of less than 0.05 will be retained. 

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_ORA.png" width="450">

### (6) Pathway analysis
Click Pathway. On the left, choose the group to analyze. All genes in the selected group will be included in the analysis and compared to the selected gene sets from the database using gene set enrichment analysis (GSEA) or other methods. Pathways with an adjusted p-value (FDR) below the chosen threshold will be retained. 

<img src="https://github.com/TK-CamBaz/Transcriptomic-data-analysis-DEMO/blob/main/Tutorial/Figure/iDEP_GSEA.png" width="450">

**Other analysis methods (Genome, Bicluster, and Network) are also available but are not listed here. Feel free to explore them if interested. 
