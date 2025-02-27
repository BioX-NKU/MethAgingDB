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
output <- 'output/'

counts_meta <- read_and_merge_files(gse_id, base_path)

counts_meta$all_meta$Group1 <-  factor(ifelse(counts_meta$all_meta$Age_group == "<20", "<20", "other"))
counts_meta$all_meta$Group2 <-  factor(ifelse(counts_meta$all_meta$Age_group == "20_40", "20_40", "other"))
counts_meta$all_meta$Group3 <-  factor(ifelse(counts_meta$all_meta$Age_group == "40_60", "40_60", "other"))
counts_meta$all_meta$Group4 <-  factor(ifelse(counts_meta$all_meta$Age_group == "60_80", "60_80", "other"))
counts_meta$all_meta$Group5 <-  factor(ifelse(counts_meta$all_meta$Age_group == ">80", ">80", "other"))
my_dict <- list(Group1 = "<20", Group2 = "20_40", Group3 = "40_60", Group4 = "60_80", Group5 = ">80")

group_names <- paste("Group", 1:5, sep="")

combined_DMP <- data.frame()
for (group_name in group_names){
        if (length(unique(counts_meta$all_meta[[group_name]])) <= 1) {
          # If the group has only one unique value, skip the current iteration
          next
        }

        # Use tryCatch to capture potential errors
        myDMP <- tryCatch({
            champ.DMP(beta = counts_meta$all_counts, pheno=counts_meta$all_meta[[group_name]])
        }, error = function(e) {
            # If an error occurs, print the error message and return NULL
            print(paste("Error in processing", group_name, ":", e$message))
            return(NULL)
        })

        # Check if myDMP is NULL, and if so, skip the current iteration
        if (is.null(myDMP)) {
            next
        }

        myDMP[[1]]$DMP_type <- my_dict[[group_name]]

        # Use dplyr's rename_with to update column names
        myDMP[[1]] <- myDMP[[1]] %>%
          rename_with(~ "Group_AVG", .cols = sprintf("X%s_AVG", my_dict[[group_name]]))
        myDMP[[1]]$value <- myDMP[[1]]$Group_AVG - myDMP[[1]]$other_AVG

        combined_DMP <- rbind(combined_DMP, myDMP[[1]])
    }
    write.csv(combined_DMP, sprintf("%s/DMP/%s.csv", base_path, gse_id))

