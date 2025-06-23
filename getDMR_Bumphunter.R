library(dplyr)
library(data.table)
library(ChAMP)
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPseeker)
library(bumphunter)
library(org.Hs.eg.db)

read_and_merge_files <- function(gse_id, base_path) {
  
  
    count_file_path <- file.path(base_path,  paste0(gse_id, "_count.csv"))
    meta_file_path <- file.path(base_path,  paste0(gse_id, "_meta.csv"))
    
    counts <- fread(count_file_path)
    meta <- fread(meta_file_path)
    rownames(counts) <- counts$pos
    counts$pos <- NULL
    rownames(meta) <- meta$GSM_number
    meta$GSM_number <- NULL
    meta$Age_group <- gsub("-", "_", meta$Age_group)
    
    return(list(all_counts = counts, all_meta = meta))
}

base_path <- '/data/'
gse_id <- 'GSE17448'
output_dir <- 'output/'

counts_meta <- read_and_merge_files(gse_id, base_path)

counts_meta$all_meta$Group1 <-  factor(ifelse(counts_meta$all_meta$Age_group == "<20", "<20", "other"))
counts_meta$all_meta$Group2 <-  factor(ifelse(counts_meta$all_meta$Age_group == "20_40", "20_40", "other"))
counts_meta$all_meta$Group3 <-  factor(ifelse(counts_meta$all_meta$Age_group == "40_60", "40_60", "other"))
counts_meta$all_meta$Group4 <-  factor(ifelse(counts_meta$all_meta$Age_group == "60_80", "60_80", "other"))
counts_meta$all_meta$Group5 <-  factor(ifelse(counts_meta$all_meta$Age_group == ">80", ">80", "other"))
my_dict <- list(Group1 = "<20", Group2 = "20_40", Group3 = "40_60", Group4 = "60_80", Group5 = ">80")

group_names <- paste("Group", 1:5, sep="")

combined_DMR <- data.frame()

for (group_name in group_names){
    
    # Check the number of unique values in the corresponding group_name column in the meta data
    if (length(unique(counts_meta$all_meta[[group_name]])) <= 1) {
      # If the group has only one unique value, skip the current iteration
      next
    }
    
    # Use tryCatch to capture any potential errors
    myDMR <- tryCatch({
        champ.DMR(beta = data.matrix(counts_meta$all_counts), pheno=counts_meta$all_meta[[group_name]], minProbes=3, adjPvalDmr=1)
    }, error = function(e) {
        # If an error occurs, print the error message and return NULL
        print(paste("Error in processing", group_name, ":", e$message))
        return(NULL)
    })
    
    # Check if myDMR is NULL, and if so, skip the current iteration
    if (is.null(myDMR)) {
        next
    }
    
    myDMR[[1]]$DMR_type <- my_dict[[group_name]]
    
    dmrs_gr <- GRanges(
      seqnames = Rle(myDMR$BumphunterDMR$seqnames),
      ranges = IRanges(start = myDMR$BumphunterDMR$start, end = myDMR$BumphunterDMR$end)
    )

    annotation <- annotatePeak(dmrs_gr, tssRegion = c(-3000, 3000), TxDb = txdb)
    annotation <- as.data.frame(annotation)
    
    # Assuming geneIds is the vector of gene IDs obtained from ChIPseeker annotation
    geneIds <- as.character(annotation$geneId)

    # Retrieve gene names
    geneNames <- mapIds(org.Hs.eg.db, keys = geneIds, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    annotation$geneName <- geneNames[match(annotation$geneId, names(geneNames))]
    
    rownames(annotation) <- paste(annotation$seqnames, annotation$start, annotation$end, sep = '_')
    rownames(myDMR[[1]]) <- paste(myDMR[[1]]$seqnames, myDMR[[1]]$start, myDMR[[1]]$end, sep = '_')
    
    myDMR[[1]] <- cbind(annotation, myDMR[[1]])
    combined_results <- rbind(combined_results, myDMR[[1]])

}
    write.csv(combined_DMR, sprintf("%s/DMR/%s.csv", output_dir, gse_id))

