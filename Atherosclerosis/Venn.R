
if (!require("VennDiagram")) install.packages("VennDiagram")
if (!require("Cairo")) install.packages("Cairo")
install.packages("Cairo")
if (!require("extrafont")) install.packages("extrafont")
library(VennDiagram)
library(Cairo)
library(extrafont)
library(grid)  




setwd("/home/xiaoning/Documents/projects/MD-sun/PFOS-CKM/data/GSE120521")


deg <- read.table("PTC_vs_ATC_DEG_with_symbol.txt", header = TRUE, sep = "\t")
wgcna <- read.table("GSE120521_mapped_genes_stabilityofregion.txt", header = TRUE, sep = "\t")
pfos_ckm_genes <- read.table("gene_intersection_483.txt", header = TRUE, sep = "\t")


deg_genes <- deg$SYMBOL
wgcna_genes <- wgcna$SYMBOL
pfos_ckm_genes <- pfos_ckm_genes$PPARG

wgcna_genes <- wgcna_genes[!is.na(wgcna_genes) & wgcna_genes != ""]
deg_genes <- deg_genes[!is.na(deg_genes) & deg_genes != ""]
pfos_ckm_genes <- pfos_ckm_genes[!is.na(pfos_ckm_genes) & pfos_ckm_genes != ""]


head(wgcna_genes)
head(deg_genes)
head(pfos_ckm_genes)


venn.plot <- venn.diagram(
  x = list(
    DEG = deg_genes,
    WGCNA = wgcna_genes,
    PFOS_CKM = pfos_ckm_genes
  ),
  filename = NULL,
  col = "transparent",
  fill = c("#9BB7D4", "#B8D8B8", "#E5C1B3"), 
  alpha = 0.5,
  cex = 1.2,  
  fontface = "bold",
  fontfamily = "Times New Roman",  
  cat.cex = 1.2,  
  cat.fontfamily = "Times New Roman",  
  cat.fontface = "bold",
  cat.col = c("#2A5CAA", "#5A8F5A", "#C76A49"),
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.09, 0.09, 0.08),
  margin = 0.15
)

CairoPNG(
  filename = "Optimized_Venn_Diagram_stabilityofregion_483.png", 
  width = 1000, 
  height = 1000, 
  res = 200
)
grid.draw(venn.plot)
dev.off()


intersection_genes <- Reduce(intersect, list(deg_genes, wgcna_genes, pfos_ckm_genes))
length(intersection_genes) 
head(intersection_genes) 
write.table(intersection_genes, "intersection_genes_stabilityofregion_483.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("✅ Save Venn png：Optimized_Venn_Diagram.png\n")



