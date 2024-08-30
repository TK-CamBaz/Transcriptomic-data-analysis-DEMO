詳細分析流程
@@ 於下方提到的軟體或擴充套件如遇安裝問題，請自行解決。

1. 準備資料
在此專案中，須先準備3份檔案，分別為RNA-seq資料、全基因組序列以及全基因組註解。資料來源為NCBI (accession: PRJNA946758)以及Figshare (https://figshare.com/articles/dataset/Chromosome-level_genome_assembly_of_the_two-spotted_spider_mite_i_Tetranychus_urticae_i_/25241794/1)。RNA-seq資料包含12個sra檔，包含各樣本的轉錄體資料，其試驗設計為3個處理組:不接觸藥劑的感性品系、持續接觸藥劑的抗性品系以及不接觸藥劑的抗性品系，每個處理有4個重複，其中感性品系為倫敦對照品系(London reference strain)，抗性品系為FP9品系(對亞醌蟎acequinocyl及必芬蟎bifenazate具高度抗性)，接觸藥劑為亞醌蟎。全基因組序列為二點葉蟎的基因體序列，組裝至染色體等級，為fasta檔。全基因組註解為二點葉蟎之編碼序列(CDS)在基因體上的位置，為gff檔。sra檔可用sra-tookit中的prefetch進行下載：

先將欲下載的sra檔名寫入文件檔(如sra-list.txt)，接著用指令
"""
prefetch --option-file sra_list.txt &
"""
# "&" 代表使這段指令於背景運行

而全基因組序列及全基因組註解則直接至該網址下載即可。

2. 資料前處理
