library(dplyr)
library(data.table)
library(ChAMP)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(bumphunter)
library(org.Hs.eg.db)


base_path <- './'
output_dir <- './output'
  
gse_ids <- 'GSE17448'
    
count_file_path <- file.path(base_path, gse_id, paste0(gse_id, "_count_matrix.csv"))
meta_file_path <- file.path(base_path,gse_id, paste0(gse_id, "_meta.csv"))
counts <- fread(count_file_path, sep = '\t')
cpg <- counts[[1]]
counts <- counts[,-1]
rownames(counts) <- cpg
meta <- read.csv(meta_file_path)

counts_meta <- list(all_counts = counts, all_meta = meta)

mat <- as.matrix(counts_meta$all_counts)
rownames(mat) <- rownames(counts_meta$all_counts)
colnames(mat) <- NULL  # 删除列名


myDMP <-  champ.DMP(beta = mat, pheno=as.numeric(counts_meta$all_meta$Age),adjPVal = 0.05)

        

        