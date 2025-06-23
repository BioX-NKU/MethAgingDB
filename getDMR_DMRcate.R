library(dplyr)
library(data.table)
library(ChAMP)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(bumphunter)
library(org.Hs.eg.db)

mychampDMR <- function (beta = myNorm, pheno = myLoad$pd$Sample_Group, compare.group = NULL, 
    arraytype = "450K", method = "Bumphunter", minProbes = 7, 
    adjPvalDmr = 0.05, cores = 3, maxGap = 300, cutoff = NULL, 
    pickCutoff = TRUE, smooth = TRUE, smoothFunction = loessByCluster, 
    useWeights = FALSE, permutations = NULL, B = 250, nullMethod = "bootstrap", 
    meanLassoRadius = 375, minDmrSep = 1000, minDmrSize = 50, 
    adjPvalProbe = 0.05, Rplot = T, PDFplot = T, resultsDir = "./CHAMP_ProbeLasso/", 
    rmSNPCH = T, fdr = 0.05, dist = 2, mafcut = 0.05, lambda = 1000, 
    C = 2) 
{
    message("[===========================]")
    message("[<<<<< ChAMP.DMR START >>>>>]")
    message("-----------------------------")
    message("!!! important !!! We just upgrate champ.DMR() function, since now champ.DMP() could works on multiple phenotypes, but ProbeLasso can only works on one DMP result, so if your pheno parameter contains more than 2 phenotypes, and you want to use ProbeLasso function, you MUST specify compare.group=c(\"A\",\"B\"). Bumphunter and DMRcate should not be influenced.")
    message("\n[ Section 1:  Check Input Pheno Start ]\n")
    if (length(which(is.na(beta))) > 0) 
        message(length(which(is.na(beta))), " NA are detected in your beta Data Set, which may cause fail or uncorrect of SVD analysis. You may want to impute NA with champ.impute() function first.")
    if (!class(pheno) %in% c("character", "factor", "numeric")) 
        stop("pheno parameter must be a category vector, could be character, factor or numeric (numeric is not recommended).")
    message("  You pheno is ", class(pheno), " type.")
    message("    Your pheno information contains following groups. >>")
    sapply(unique(pheno), function(x) message("    <", x, ">:", 
        sum(pheno == x), " samples."))
   if (method == "DMRcate") {
        message(cores, " cores will be used to do parallel DMRcate computing.")
        message("<< Find DMR with DMRcate Method >>")
        myMs <- logit2(beta)
        if (rmSNPCH) 
            myMs <- rmSNPandCH(myMs, dist = dist, mafcut = mafcut)
        design <- model.matrix(~pheno)
        if (arraytype == "450K") {
            myannotation <- cpg.annotate(datatype = "array", 
                fdr = fdr, myMs, design = design, coef = ncol(design), 
                analysis.type = "differential", annotation = c(array = "IlluminaHumanMethylation450k", 
                  annotation = "ilmn12.hg19"), what = "M")
        }
        else {
            myannotation <- cpg.annotate(datatype = "array", 
                fdr = fdr, myMs, design = design, coef = ncol(design), 
                analysis.type = "differential", annotation = c(array = "IlluminaHumanMethylationEPIC", 
                  annotation = "ilm10b4.hg19"), what = "M")
        }
#         M <- do.call("cbind", lapply(myannotation, as.data.frame))
#         colnames(M) <- names(myannotation)
        dmrcoutput <- dmrcate(myannotation, min.cpgs = minProbes, 
            lambda = lambda, C = C)
        data(dmrcatedata)
        DMR <- as.data.frame(extractRanges(dmrcoutput, genome = "hg19"))
        message("DMRcate detected ", nrow(DMR), " DMRs with mafcut as= ", 
            adjPvalDmr, ".")
        if (nrow(DMR) == 0) 
            stop("No DMR detected.")
        if (nrow(DMR) != 0) {
            rownames(DMR) <- paste("DMR", 1:nrow(DMR), sep = "_")
            OutputDMR <- list(DMRcateDMR = DMR)
        }
        else {
            OutputDMR <- NULL
        }
    }
    else {
        stop("Please assign correct DMR method: 'Bumphunter' or 'ProbeLasso'")
    }
    message("\n[ Section 2:  Run DMR Algorithm Done ]\n")
    message("[<<<<<< ChAMP.DMR END >>>>>>]")
    message("[===========================]")
    message("[You may want to process DMR.GUI() or champ.GSEA() next.]\n")
    return(OutputDMR)
}

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
meta$Age_group <- gsub("-", "_", meta$Age_group)
meta$Age_group <- gsub("<", "lt", meta$Age_group)
meta$Age_group <- gsub(">", "mt", meta$Age_group)

counts_meta <- list(all_counts = counts, all_meta = meta)

counts_meta$all_meta$Group1 <-  factor(ifelse(counts_meta$all_meta$Age_group == "lt20", "lt20", "other"))
counts_meta$all_meta$Group2 <-  factor(ifelse(counts_meta$all_meta$Age_group == "20_40", "20_40", "other"))
counts_meta$all_meta$Group3 <-  factor(ifelse(counts_meta$all_meta$Age_group == "40_60", "40_60", "other"))
counts_meta$all_meta$Group4 <-  factor(ifelse(counts_meta$all_meta$Age_group == "60_80", "60_80", "other"))
counts_meta$all_meta$Group5 <-  factor(ifelse(counts_meta$all_meta$Age_group == "mt80", "mt80", "other"))
my_dict <- list(Group1 = "lt20", Group2 = "20_40", Group3 = "40_60", Group4 = "60_80", Group5 = "mt80")
group_names <- paste("Group", 1:5, sep="")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

combined_DMR <- data.frame()

for (group_name in group_names){

    # 检查meta中对应group_name列的唯一值数量
    if (length(unique(counts_meta$all_meta[[group_name]])) <= 1) {
        # 如果该组只有一个唯一值，则跳过当前循环迭代
        next
    }

    # 使用tryCatch捕获可能发生的错误
    myDMR <- tryCatch({
        mychampDMR(beta = data.matrix(counts_meta$all_counts), pheno=counts_meta$all_meta[[group_name]],minProbes=3,method='DMRcate')
    }, error = function(e) {
        # 如果发生错误，打印错误消息并返回NULL
        print(paste("Error in processing", group_name, ":", e$message))
        return(NULL)
    })
        # 检查myDMR是否为NULL，如果是，则跳过当前迭代
    if (is.null(myDMR)) {
        next
    }
    myDMR[[1]]$DMR_type <- my_dict[[group_name]]

    dmrs_gr <- GRanges(
        seqnames = Rle(myDMR$DMRcateDMR$seqnames),
        ranges = IRanges(start = myDMR$DMRcateDMR$start, end = myDMR$DMRcateDMR$end)
    )

    annotation <- annotatePeak(dmrs_gr, tssRegion = c(-3000, 3000), TxDb = txdb)
    annotation <- as.data.frame(annotation)
    # 假设 geneIds 是你从 ChIPseeker 注释得到的基因ID向量
    geneIds <- as.character(annotation$geneId)

    # 获取基因名称
    geneNames <- mapIds(org.Hs.eg.db, keys = geneIds, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    annotation$geneName <- geneNames[match(annotation$geneId, names(geneNames))]
    rownames(annotation) <- paste(annotation$seqnames,annotation$start,annotation$end,sep = '_')
    rownames(myDMR[[1]]) <- paste(myDMR[[1]]$seqnames,myDMR[[1]]$start,myDMR[[1]]$end,sep = '_')
    myDMR[[1]] <- cbind(annotation,myDMR[[1]])
    myDMR[[1]]$DMR_name <- rownames(myDMR[[1]])
    rownames(myDMR[[1]]) <- NULL
    combined_DMR <- rbind(combined_DMR, myDMR[[1]])
    }

write.csv(combined_DMR, sprintf("%s/DMRcate/%s_%s.csv", output_dir, tissue,gse_id))


}}