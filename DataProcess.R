# Load annotation file
require(data.table)
emapper.anno <- fread('Partial.emapper.annotations.tsv', sep = '\t')

anno.id <- c()
for (i in 1:nrow(emapper.anno)){
  temp.anno.id <- strsplit(emapper.anno$seed_ortholog[i], '.', fixed = T)
  temp.anno.id <- temp.anno.id[[1]][2]
  anno.id <- c(anno.id, temp.anno.id)
}

emapper.anno$gene_id <- anno.id
colnames(emapper.anno)[1] <- 'query'
View(emapper.anno)

# Load read counts file
require(data.table)
read.count <- fread('Partial.Read.Counts.csv', sep = ',')
colnames(read.count)[1] <- 'query'
read.count$query <- gsub('TU', 'model', read.count$query)

# Merge two files
anno.read.merge <- merge(emapper.anno[, c('query', 'gene_id')], read.count, by = 'query')
read.count.for.iDEP <- anno.read.merge[, -1]
colnames(read.count.for.iDEP) <- c('gene_id', 'SRR23919708', 'SRR23919707', 'SRR23919704',
                                   'SRR23919703', 'SRR23919698', 'SRR23919697',
                                   'SRR23919706', 'SRR23919705', 'SRR23919702',
                                   'SRR23919701', 'SRR23919700', 'SRR23919699')

# Save file for follow-up analysis
write.csv(read.count.for.iDEP, 'read.count.for.iDEP.csv', quote = F, row.names = F)
