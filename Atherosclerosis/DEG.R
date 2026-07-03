setwd("/home/xiaoning/Documents/projects/MD-sun/PFOS-CKM/data/GSE120521")
library(DESeq2)
dataMatrix<-read.table("GSE120521_raw_counts_GRCh38.p13_NCBI.tsv",header = T)
dim(dataMatrix)
head(dataMatrix)
rownames(dataMatrix)<-dataMatrix[,"GeneID"]
dataMatrix<-dataMatrix[,2:dim(dataMatrix)[2]]
head(dataMatrix)
coldata <- read.table("GSE120521_deg_sampleinfo.txt",header = T)
head(coldata)
colnames(coldata)
dds <- DESeqDataSetFromMatrix(countData = dataMatrix, colData = coldata, design = ~ stabilityofregion)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
dds.res <- results(dds,alpha = 0.05)
summary(dds.res)
sum(dds.res$padj<0.05, na.rm=TRUE)
dds.resOrdered <- dds.res[order(dds.res$padj),]
write.table(as.data.frame(dds.resOrdered), file="PTC_vs_ATC_results.txt",quote = F,sep = "\t")
sum( (dds.res$padj<0.05 & abs(dds.res$log2FoldChange)>2),na.rm=TRUE)
save(dds.res,file="8.5.dds.res.Rdata") 

library(AnnotationDbi)
library(org.Hs.eg.db)
res <- read.table(
  "PTC_vs_ATC_results.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE
)
head(rownames(res))
# 3151 107984950 9429 8932
ids <- rownames(res)
mapped_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = ids,
  columns  = "SYMBOL",
  keytype  = "ENTREZID"
)
mapped_genes <- mapped_genes[!duplicated(mapped_genes$ENTREZID), ]
mapped_genes <- mapped_genes[match(ids, mapped_genes$ENTREZID), ]
head(mapped_genes)
res$SYMBOL <- mapped_genes$SYMBOL
res <- res[, c("SYMBOL", colnames(res)[1:(ncol(res)-1)])]
write.table(
  res,
  file = "PTC_vs_ATC_results_with_SYMBOL.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)



dds.res.filtered <- dds.res[which(dds.res$padj<0.05 & abs(dds.res$log2FoldChange) >2),]
write.table(dds.res.filtered, file="PTC_vs_ATC_DEG.txt",quote = F,sep = "\t")


if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")
}
library(org.Hs.eg.db)
library(AnnotationDbi)


data <- read.table("PTC_vs_ATC_DEG.txt", header = TRUE, sep = "\t", check.names = FALSE)


first_id <- rownames(data)[1] 

if (grepl("^ENSG", first_id, ignore.case = TRUE)) {
  id_type <- "ENTREZID"
  message("Detected Ensembl Gene IDs.")
} else if (grepl("^[0-9]+$", first_id)) {
  id_type <- "ENTREZID"
  message("Detected Entrez Gene IDs.")
} else {
  stop("❌ Could not detect ID type automatically. Please check the first column.")
}


ids <- rownames(data)
mapped_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ids,     # Your vector of Entrez IDs
  columns = "SYMBOL",   # The information you want to retrieve
  keytype = "ENTREZID"  # Specify that the keys are Entrez IDs
)


# Remove duplicates if multiple mappings
mapped_genes <- mapped_genes[!duplicated(mapped_genes$ENTREZID), ]
# ===== Step 4: Merge conversion results =====

colnames(mapped_genes) <- c("ENTREZID", "SYMBOL")

data$ENTREZID <- rownames(data)  
#data$ENTREZID <- data[[1]]


data$ENTREZID <- as.character(data$ENTREZID)
mapped_genes$ENTREZID <- as.character(mapped_genes$ENTREZID)


data_merged <- merge(data, mapped_genes, by = "ENTREZID", all.x = TRUE)


head(data_merged)



data_merged <- data_merged[, c("SYMBOL", colnames(data))]


write.table(data_merged, "PTC_vs_ATC_DEG_with_symbol.txt", sep = "\t", quote = FALSE, row.names = FALSE)

message("✅ Conversion complete! Output saved as: PTC_vs_ATC_DEG_with_symbol.txt")
