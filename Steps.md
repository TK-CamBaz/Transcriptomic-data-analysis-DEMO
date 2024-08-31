詳細分析流程
@@ 於下方提到的軟體或擴充套件如遇安裝問題，請自行解決。
@@ 如對指令內的引數(arguments)有疑問，也請自行查詢。

1. 準備資料
在此專案中，我準備3份檔案，分別為RNA-seq資料、全基因組序列以及全基因組註解。資料來源為NCBI (accession: PRJNA946758)以及Figshare (https://figshare.com/articles/dataset/Chromosome-level_genome_assembly_of_the_two-spotted_spider_mite_i_Tetranychus_urticae_i_/25241794/1)。RNA-seq資料包含12個sra檔，包含各樣本的轉錄體資料，定序平台為Illumina NovaSeq 6000，layout為雙端定序。試驗設計為3個處理組:不接觸藥劑的感性二點葉蟎品系、持續接觸藥劑的抗性二點葉蟎品系以及不接觸藥劑的抗性二點葉蟎品系，每個處理有4個重複，其中感性品系為倫敦對照品系(London reference strain)，抗性品系為FP9品系(對亞醌蟎acequinocyl及必芬蟎bifenazate具高度抗性)，接觸藥劑為亞醌蟎。全基因組序列為二點葉蟎的基因體序列，組裝至染色體等級，為fasta檔。全基因組註解為二點葉蟎之編碼序列(CDS)在基因體上的位置，為gff檔。sra檔可用sra-tookit中的prefetch進行下載：
以SRX19731018為例：
"""
prefetch SRX19731018
"""
而全基因組序列及全基因組註解則直接至該網址下載即可。

2. 資料前處理
(1)將sra檔內的fastq檔取出並處理乾淨
首先將sra檔解壓縮並使用sra-tookit中的parallel-fastq-dump提取裡面的fastq檔：
以SRX19731018.sra為例：
"""
parallel-fastq-dump -s SRX19731018.sra -t 8 --split-3
"""
樣本為雙端定序的情況下會得到兩個fastq檔，兩個都要保留。

接著利用AfterQC清洗序列，主要目的為移除轉接頭(adapter)以及過短的序列，而使用AfterQC的好處在於可自動判斷清洗條件並輸出清洗序列前後的報表。
以SRR23919699為例：
"""
cd ~/AfterQC-0.9.7
pypy after.py -1 SRR23919699_1.fastq -2 SRR23919699_2.fastq
"""
輸出資料包括三個資料夾：good、bad及QC，good中的fq檔即為清洗後的乾淨序列，bad中的fq檔為被移除的序列，QC內的html檔為各序列清理前後的報表，可檢查序列清理前後的差別。

(2)建立全基因組索引
目的為加快後續序列比對的速度，使用HISAT2來執行。
以tu.genome.ipm_v2.fasta為例：
"""
hisat2-build -p 8 tu.genome.ipm_v2.fasta index/tu.genome_index
"""
輸出資料為名為tu.genome_index的資料夾，裡面有8個.h2檔。

(3)將序列比對到參考基因組上
使用HISAT2來執行。
以SRR23919699_1.good.fq及SRR23919699_2.good.fq為例：
hisat2 -t -p 2 -x index/tu.genome_index -1 SRR23919699_1.good.fq -2 SRR23919699_2.good.fq -S SRR23919699.sam
輸出資料為sam檔。

接著將sam檔轉換成bam檔，此一步驟為格式轉換，在不改寫內容的前提下大幅縮減檔案大小。
使用samtools來執行。
以SRR23919699.sam為例：
samtools sort -@2 -m 200M -o SRR23919699.bam SRR23919699.sam
輸出資料為bam檔。

接著可用qualimap檢查比對結果
"""
cd ~/qualimap_v2.3
./qualimap
"""
待視窗開啟後，選取bam檔進行分析，完成後會顯示報表，亦可輸出成pdf或html。

另外還可用bam檔生成bai檔，此檔案可與bam檔一起輸入IGV(Integrative Genomics Viewer)中將比對結果視覺化呈現：
samtools index SRR23919699.bam SRR23919699.bai
接著打開IGV
"""
igv
"""
選取bam檔和全基因組序列，注意bam檔與bai檔須放在同一路徑下。

(4)計算各樣本的轉錄本在序列比對後的數量
該步驟需要bam檔與gtf檔，先對全基因組註解的gff檔使用gffread進行轉換：
以tu_evm_out.gff3為例：
gffread tu_evm_out.gff3 -T -o tu_evm_out.gtf

接著執行HTSeq中的htseq-count：
以SRR23919699.bam為例：
"""
htseq-count -q -f bam -s no -i gene_id SRR23919699.bam tu_evm_out.gtf > SRR23919699.count
"""
重複上述步驟，得到所有樣本的count檔，利用excel或其他軟體(R或python等)將檔案合併。

3. 序列代號轉換成基因代號
**若gff檔中的ID已經是有效的基因代號，則可省略此步驟。

在合併count檔後獲得的文件檔中，第一欄為轉錄本的序列代號，並非有效的基因代號(Gene ID)，因此需要先利用線上資料庫將序列代號轉換為基因代號。首先用全基因組序列和全基因組註解(這裡要用gtf)將各轉錄本的序列取出：
"""
gffread -w transcripts.fa -g tu.genome.ipm_v2.fasta tu_evm_out.gtf
"""

接著使用eggNOG-mapper(http://eggnog-mapper.embl.de/)，選擇CDS後將轉錄本序列檔案(transcripts.fa)上傳，輸入欲收到檔案的email網址後點選提交(Submit)，待比對完成後下載csv檔(選取csv，下載後的檔案副檔名是tsv，但不影響檔案內容)。