# Analysis Workflow 
## 1. Data Preparation
In this project, I prepared three files: transcriptome data, whole genome sequences, and whole genome annotations. The data sources are from NCBI (accession: PRJNA946758) and Figshare (https://figshare.com/articles/dataset/Chromosome-level_genome_assembly_of_the_two-spotted_spider_mite_i_Tetranychus_urticae_i_/25241794/1). The transcriptome data includes 12 sra files, containing the transcriptome data of each sample. The sequencing platform is Illumina NovaSeq 6000, and the layout is paired-end sequencing. The experimental design includes three treatment groups: the susceptible two-spotted spider mite strain without exposure to pesticides, the resistant two-spotted spider mite strain with continuous pesticide exposure, and the resistant two-spotted spider mite strain without pesticide exposure, with four replicates for each treatment. The susceptible strain is the London reference strain, and the resistant strain is the FP9 strain (highly resistant to acequinocyl and bifenazate). The pesticide used is acequinocyl. The whole genome sequence is the genome assembly of the two-spotted spider mite, assembled to the chromosome level, in fasta format. The whole genome annotation is the location of the coding sequences (CDS) of the two-spotted spider mite on the genome, in gff format. The sra files can be downloaded using the prefetch command from the sra-toolkit: For example, with SRX19731018: 
```
prefetch SRX19731018
```
The whole genome sequence and annotations can be downloaded directly from the provided link.

## 2. Data Preprocessing
### (1) Extract fastq files from sra files and perform quality control
First, decompress the sra files and use the parallel-fastq-dump tool from the sra-toolkit to extract the fastq files:
For example, with SRX19731018.sra:
```
parallel-fastq-dump -s SRX19731018.sra -t 8 --split-3
```
In the case of paired-end sequencing, two fastq files will be generated, both of which should be retained.

Next, use AfterQC to clean the sequences, mainly to remove adapters and sequences that are too short. AfterQC automatically determines the cleaning conditions and outputs reports comparing the sequences before and after cleaning.
For example, with SRR23919699:
```
cd ~/AfterQC-0.9.7
pypy after.py -1 SRR23919699_1.fastq -2 SRR23919699_2.fastq
```
The output includes three directories: good, bad, and QC. The fq files in the good folder are the cleaned sequences, the fq files in the bad folder are the removed sequences, and the html files in the QC folder are the reports that compare the sequences before and after cleaning.

### (2) Build a whole genome index
The purpose of this step is to speed up the subsequent sequence alignment process, which is done using HISAT2.
For example, with tu.genome.ipm_v2.fasta:
```
hisat2-build -p 8 tu.genome.ipm_v2.fasta index/tu.genome_index
```
The output is a directory named tu.genome_index, containing eight .h2 files.

### (3) Align sequences to the reference genome
Use HISAT2 for this step.
For example, with SRR23919699_1.good.fq and SRR23919699_2.good.fq:
```
hisat2 -t -p 2 -x index/tu.genome_index -1 SRR23919699_1.good.fq -2 SRR23919699_2.good.fq -S SRR23919699.sam
```
The output is a sam file.

Next, convert the sam file to a bam file. This step is a format conversion that significantly reduces the file size without altering its content.
Use samtools for this step.
For example, with SRR23919699.sam:
```
samtools sort -@2 -m 200M -o SRR23919699.bam SRR23919699.sam
```
The output is a bam file.

Then, use qualimap to check the alignment results: 
```
cd ~/qualimap_v2.3
./qualimap
```
After the window opens, select the bam file for analysis. Once complete, a report will be displayed, which can also be exported as a pdf or html.

You can also generate a bai file from the bam file, which can be loaded into IGV (Integrative Genomics Viewer) along with the bam file for visualizing the alignment results: 
```
samtools index SRR23919699.bam SRR23919699.bai
```
Then open IGV:
```
igv
```
Select the bam file and the whole genome sequence. Make sure the bam file and bai file are in the same directory.

### (4) Quantify the transcripts in each sample after sequence alignment
This step requires the bam file and a gtf file. First, convert the whole genome annotation gff file to a gtf file using gffread:
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

## 3. Analyze the count matrix using iDEP
iDEP is a web tool designed for RNA-seq data analysis, offering common analysis methods such as data exploration, differential gene expression analysis, over-representation analysis, and gene set enrichment analysis. Besides uploading the count matrix, you can also prepare and upload an experimental design file. The format can be checked by clicking Info. In this experiment, the material used is the two-spotted spider mite, so make sure to select the correct species in the Species section. 
### (1) Preprocessing
Click Pre-Process, and each button allows you to view the distribution of Read counts across samples. Selecting rlog for expression value transformation helps to improve the correlation among samples under the same treatment. 
### (2) 聚Clustering analysis
Click Clustering. Two options are available: hierarchical clustering and k-means clustering. The clustering results are presented as a heatmap showing the gene expression profiles of each sample. You can select regions of interest on the heatmap to generate a sub-heatmap, which allows further investigation of the expression profiles. The sub-heatmap also lists gene names for easy tracking and recording.
### (3) Principal component analysis
Click PCA to observe the similarity between samples. Other dimensionality reduction methods are also available.  
### (4) Differential gene expression analysis
Click DEG1, select the comparison groups, and click Submit. The results are displayed as bar plots and tables showing the number of upregulated and downregulated differentially expressed genes in each comparison. The results can be downloaded from the Result & data section. 
### (5) Differential gene expression visualization and over-representation analysis
Click DEG2. You can observe the distribution of differentially expressed genes across the genome using heatmaps, volcano plots, MA plots, and scatter plots. Select Enrichment to perform over-representation analysis. On the left, choose the group to analyze. The differentially expressed genes in this group will be treated as a gene list and compared to the selected gene sets from the database using a hypergeometric test. Pathways with an adjusted p-value (FDR) of less than 0.05 will be retained. 
### (6) Pathway analysis
Click Pathway. On the left, choose the group to analyze. All genes in the selected group will be included in the analysis and compared to the selected gene sets from the database using gene set enrichment analysis (GSEA) or other methods. Pathways with an adjusted p-value (FDR) below the chosen threshold will be retained. 

**Other analysis methods (Genome, Bicluster, and Network) are also available but are not listed here. Feel free to explore them if interested. 