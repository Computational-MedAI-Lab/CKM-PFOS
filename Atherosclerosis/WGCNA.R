#BiocManager::install("WGCNA")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("BiocManager")
# BiocManager::install("preprocessCore", force = TRUE)
library(preprocessCore)
# BiocManager::install("WGCNA")
library(WGCNA)
#install.packages("WGCNA")

setwd("/home/xiaoning/Documents/projects/MD-sun/PFOS-CKM/data/GSE120521")
WGCNA_matrix <- read.table("GSE120521_norm_counts_FPKM_GRCh38.p13_NCBI.tsv", 
                           sep = "\t", 
                           header = TRUE, 
                           row.names = 1)

sampleInfor <- read.table("GSE120521_sampleinfo.txt", sep = "\t", fill = TRUE, header = TRUE, row.names = 1)
head(sampleInfor)

#datExpr0 <- t(WGCNA_matrix[order(apply(WGCNA_matrix,1,var), decreasing = T)[1:2000],])  # select var proteins
datExpr0 <- t(WGCNA_matrix[order(apply(WGCNA_matrix,1,var), decreasing = T)[1:2000],]) 
gsg = goodSamplesGenes(datExpr0, verbose = 3)   
gsg$allOK


datExpr <- datExpr0
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf("GSE120521.pdf", width=36, height=12)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()



nGenes = ncol(datExpr)  
nSamples = nrow(datExpr)  

datTraits <- sampleInfor[colnames(WGCNA_matrix),]


powers = c(c(1:10), seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
power = sft$powerEstimate
power 


# ===============================
powers <- c(1:30)

options(mc.cores = 1)
disableWGCNAThreads()

sft <- pickSoftThreshold(datExpr, powerVector = powers, blockSize = 5000, verbose = 5)

write.csv(sft$fitIndices, "SoftThresholdingResults.csv", row.names = FALSE)

pdf("SoftThresholdPlots_20260630.pdf", width = 12, height = 6, family = "Times")  

par(mfrow = c(1, 2))
# -------------------------------
# 1) Scale independence
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.75, col = "red")
# -------------------------------
# 2) Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = 0.9, col = "red")

dev.off()
power=13  


net = blockwiseModules(datExpr, power = power,TOMType = "unsigned", minModuleSize = 20,reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "TOM",verbose = 3)
table(net$colors)  



mergedColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]]
pdf("GSE120521.pdf", width=15, height=10, family = "Times")
plotDendroAndColors(geneTree, mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(stringr::str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col) 


pdf("GSE120521.pdf", width=10, height=12)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()


design=model.matrix(~0+ datTraits$stabilityofregion)
colnames(design)=levels(datTraits$stabilityofregion)
moduleTraitCor = cor(MEs_col, design, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")  # signif表示保留几位小数
dim(textMatrix) = dim(moduleTraitCor)
sizeGrWindow(9,9)
par(mar = c(3, 8, 3, 1));
unique_labels <- unique(datTraits$stabilityofregion)
colnames(moduleTraitCor) <- unique_labels
label_map <- c(
  stable   = "Stable",
  unstable = "Unstable"
)
xLabels_plot <- label_map[colnames(moduleTraitCor)]
pdf("stabilityofregion.pdf",
    width = 7, height = 10, family = "Times") 

par(mar = c(8, 8, 4, 2) + 0.1)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = xLabels_plot,
  yLabels = names(MEs_col),
  ySymbols = names(MEs_col),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  
  cex.text  = 1.0,
  cex.lab.x = 1.3,
  xLabelsAngle = 0,
  xLabelsAdj   = 0.5, 
  
  zlim = c(-1, 1),
  main = "Module–trait relationships"
)

dev.off()



moduleColors = labels2colors(net$colors)  

table(moduleColors)

pNames = colnames(datExpr) 
# Define the list of modules of interest by color
modules_of_interest <- c("brown", "blue", "pink", "turquoise", "black")

# Create a logical vector indicating which modules are of interest
inModule <- moduleColors %in% modules_of_interest

# Extract the pNames corresponding to the modules of interest
modProbes <- pNames[inModule]

# Check the length of the selected probes
length(modProbes)

# View the extracted pNames (optional)
modProbes

write.table(modProbes, file = "modProbes_stabilityofregion.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

library(org.Hs.eg.db)


modProbes <- read.table("modProbes_stabilityofregion.txt", header = FALSE, stringsAsFactors = FALSE)
# Check the first few elements of modProbes
head(modProbes)

modProbes <- as.character(modProbes)

typeof(modProbes)  
mapped_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = modProbes,      
  columns = "SYMBOL",     
  keytype = "ENTREZID"   
)
head(mapped_genes)

missing_genes <- mapped_genes[is.na(mapped_genes$SYMBOL), ]

mapped_genes_clean <- mapped_genes[!is.na(mapped_genes$SYMBOL), ]

head(mapped_genes_clean)

write.table(mapped_genes, "GSE55296_mapped_genes_stabilityofregion.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


