library(clusterProfiler)   库(clusterProfiler)
library(org.Hs.eg.db)   库(org.Hs.eg.db)
library(ggplot2)   库(ggplot2)
library(ggnewscale)   图书馆(ggnewscale)
library(stringr)
library(scales)   库(尺度)
library(DOSE)   库(剂量)
library(enrichplot)
library(tidyverse)   库(tidyverse)
library(openxlsx)   库(openxlsx)

# Set working directory
path="C:/Users/29077/Desktop/20251203-enrichment"
setwd(path)   setwd(路径)

# Load gene list (Entrez Gene ID)
genes_raw <- read.table("gene_intersection_482.txt", header = FALSE)$V1
gene_list <- genes_raw[genes_raw != "gene_intersection_482"]
length(gene_list)   长度(gene_list)

# Convert Gene Symbol to Entrez ID (required for KEGG analysis)
gene_entrez <- mapIds(org.Hs.eg.db,gene_entrez <- mapIds(org. hs . egg .db)；
                      keys = gene_list,   键= gene_list，
                      column = "ENTREZID",
                      keytype =    keytype = "SYMBOL"，"SYMBOL",
                      multiVals = "first")
gene_entrez <- na.omit(gene_entrez)

# GO enrichment analysis
ego_bp <- enrichGO(gene          = gene_entrez,ego_bp <- enrichment go (gene = gene_entrez，
                   OrgDb         = org.Hs.eg.db,
                   keyType       = keyType = 'ENTREZID'，'ENTREZID',
                   ont           = ont = "BP"，"BP",
                   pAdjustMethod =    pAdjustMethod = "BH"，"BH",
                   pvalueCutoff  =    pvalueCutoff = 0.05，0.05,
                   qvalueCutoff  =    qvalueCutoff = 0.2；0.2,
                   readable      = TRUE)

ego_mf <- enrichGO(gene          = gene_entrez,ego_mf <- enrichment go (gene = gene_entrez，
                   OrgDb         = org.Hs.eg.db,
                   keyType       = keyType = 'ENTREZID'，'ENTREZID',
                   ont           = ont = "MF"，"MF",
                   pAdjustMethod =    pAdjustMethod = "BH"，"BH",
                   pvalueCutoff  =    pvalueCutoff = 0.05，0.05,
                   qvalueCutoff  =    qvalueCutoff = 0.2；0.2,
                   readable      = TRUE)

ego_cc <- enrichGO(gene          = gene_entrez,ego_cc <- enrichment go (gene = gene_entrez，
                   OrgDb         = org.Hs.eg.db,
                   keyType       = keyType = 'ENTREZID'，'ENTREZID',
                   ont           = ont = "CC"，"CC",
                   pAdjustMethod =    pAdjustMethod = "BH"，"BH",
                   pvalueCutoff  =    pvalueCutoff = 0.05，0.05,
                   qvalueCutoff  =    qvalueCutoff = 0.2；0.2,
                   readable      = TRUE)
print(ego_bp)
print(ego_mf)
print(ego_cc)

# Combine GO results and add ontology labels
go_bp_df <- as.data.frame(ego_bp)Go_bp_df <- as.data.frame（ego_bp）
go_bp_df$Category <- go_bp_df$Category <-“生物过程”"Biological Process"

go_mf_df <- as.data.frame(ego_mf)
go_mf_df$Category <- "Molecular Function"

go_cc_df <- as.data.frame(ego_cc)
go_cc_df$Category <- "Cellular Component"

combined_go <- rbind(go_bp_df, go_mf_df, go_cc_df)

# Save results
write.xlsx(as.data.frame(combined_go), "GO_enrichment.xlsx")
write.xlsx(as.data.frame(go_bp_df), "GO_enrichment_bp.xlsx")
write.xlsx(as.data.frame(go_mf_df), "GO_enrichment_mf.xlsx")
write.xlsx(as.data.frame(go_cc_df), "GO_enrichment_cc.xlsx")

###################################################################
library(circlize)
library(grid)
library(graphics)
library(readxl)
# Read GO enrichment results
bp <- read.xlsx("GO_enrichment_bp.xlsx")
mf <- read.xlsx("GO_enrichment_mf.xlsx")
cc <- read.xlsx("GO_enrichment_cc.xlsx")

# BP: sort by qvalue
bp_sorted <- bp[order(bp$qvalue),]
bp_top_ids <- bp_sorted$ID[1:8]
print(bp_top_ids)

# MF: sort by qvalue
mf_sorted <- mf[order(mf$qvalue), ]
mf_top_ids <- mf_sorted$ID[1:6]
print(mf_top_ids)

# CC: sort by qvalue
cc_sorted <- cc[order(cc$qvalue), ]
cc_top_ids <- cc_sorted$ID[1:6]
print(cc_top_ids)

combined_go_q <- c(cc_top_ids, mf_top_ids, bp_top_ids)
print(combined_go_q)
go_ids <- as.vector(t(combined_go_q))

##############################################################################
## Add -log10(qvalue)
# BP
bp_top_qvalue <- bp_sorted$qvalue[1:8]
print(bp_top_qvalue)
bp_top_log10qvalue <- -log10(bp_top_qvalue)
print(bp_top_log10qvalue)

# MF
mf_top_qvalue <- mf_sorted$qvalue[1:6]
print(mf_top_qvalue)
mf_top_log10qvalue <- -log10(mf_top_qvalue)
print(mf_top_log10qvalue)

# CC
cc_top_qvalue <- cc_sorted$qvalue[1:6]
print(cc_top_qvalue)
cc_top_log10qvalue <- -log10(cc_top_qvalue)
print(cc_top_log10qvalue)
combine_go_log10qvalue <- c(cc_top_log10qvalue, mf_top_log10qvalue, bp_top_log10qvalue)
print(combine_go_log10qvalue)

##################################################################
# Add number of genes and number of selected genes
# Number of background genes
bp_genes_number <- sapply(strsplit(bp_sorted$BgRatio[1:8], "/"), function(x) as.numeric(x[1]))
print(bp_genes_number)

mf_genes_number <- sapply(strsplit(mf_sorted$BgRatio[1:6], "/"), function(x) as.numeric(x[1]))
print(mf_genes_number)

cc_genes_number <- sapply(strsplit(cc_sorted$BgRatio[1:6], "/"), function(x) as.numeric(x[1]))
print(cc_genes_number)

combine_number_of_genes <- c(cc_genes_number, mf_genes_number, bp_genes_number)
print(combine_number_of_genes)

# Number of selected genes
bp_select_number <- sapply(strsplit(bp_sorted$GeneRatio[1:8], "/"), function(x) as.numeric(x[1]))
print(bp_select_number)

mf_select_number <- sapply(strsplit(mf_sorted$GeneRatio[1:6], "/"), function(x) as.numeric(x[1]))
print(mf_select_number)

cc_select_number <- sapply(strsplit(cc_sorted$GeneRatio[1:6], "/"), function(x) as.numeric(x[1]))
print(cc_select_number)
combine_number_of_select <- c(cc_select_number, mf_select_number, bp_select_number)
print(combine_number_of_select)
################################################################
# Add RichFactor

bp_richfactor <- bp_sorted$RichFactor[1:8]
print(bp_richfactor)

mf_richfactor <- mf_sorted$RichFactor[1:6]
print(mf_richfactor)

cc_richfactor <- cc_sorted$RichFactor[1:6]
print(cc_richfactor)

combine_richfactor <- c(cc_richfactor, mf_richfactor, bp_richfactor)
#################################################################
cat("go_ids length:", length(go_ids), "\n")
cat("combine_go_log10qvalue length:", length(combine_go_log10qvalue), "\n")
cat("combine_number_of_genes length:", length(combine_number_of_genes), "\n")
cat("combine_number_of_select length:", length(combine_number_of_select), "\n")
cat("combine_richfactor length:", length(combine_richfactor), "\n")

# Create data frame for plotting
go_df <- data.frame(
  Number = 1:length(go_ids),
  GO_ID = go_ids,
  log10qvalue = combine_go_log10qvalue,
  number_of_genes = combine_number_of_genes, 
  number_of_select = combine_number_of_select, 
  richfactor = combine_richfactor,
  Category = rep(c("CC", "MF", "BP"), times = c(6, 6, 8)),
  stringsAsFactors = FALSE
)

print(go_df)
###############Plotting
tiff("GO_Enrichment_Circos.tiff", width = 23000, height = 23000, res = 1600, compression = "lzw")

# Global font setting
par(family = "serif", font = 1)

#########################################################################
# Colors and scales
category_colors <- c("CC" = "#7895C1", "MF" = "#A8CBDF","BP"= "#8074C8")
scale_range <- c(0, 500) 

# Circos parameters
circos.clear()
circos.par(
  start.degree = 90, 
  clock.wise = TRUE,
  track.margin = c(0, 0),
  "track.height" = 0.1
)
sector = go_df$GO_ID
s1 = factor(sector)

circos.initialize(
  s1,
  xlim = scale_range
)
# Outer track: GO IDs
circos.track(s1, ylim = c(0,1),
             bg.col = category_colors[go_df$Category],
             panel.fun = function(x,y){
               circos.text(CELL_META$xcenter,
                           CELL_META$cell.ylim[1] + 0.7,
                           #CELL_META$sector.index,
                           labels = go_df$GO_ID[CELL_META$sector.numeric.index],
                           facing = "bending.inside",
                           niceFacing = TRUE,
                           col = "black",
                           cex = 1.8,
               )
               circos.axis(labels.cex = 1.3,
                           major.tick.length = 0.35,
                           major.at = c(0, 100, 300, 500)
               )
             })

##############################################################
# Legend: Ontology categories (bottom-left)
#circos.clear()
category_color <- c(
  "Biological Process" = "#8074C8",  
  "Molecular Function" = "#A8CBDF", 
  "Cellular Component" = "#7895C1"  
)

legend(
  "bottomleft",                      
  legend = names(category_color),   
  fill = category_color,            
  title = "ONTOLOGY",               
  border = NA,                       
  bty = "n",                         
  cex = 1.8,                         
  text.font = 1,                    
  title.font = 2
)
###############################################################
# Legend: -log10(qvalue) gradient (bottom-right)
#circos.clear()

qvalue_range <- range(go_df$log10qvalue)
qvalue_breaks <- round(seq(qvalue_range[1], qvalue_range[2], length.out = 9))
qvalue_colors <- colorRampPalette(c("#F0C284", "#E3625D"))(8) 
#qvalue_breaks <- c("(0,2]", "(2,4]", "(4,6]", "(6,8]", "(8,10]", "(10,15]", "(15,20]", ">=20")
qvalue_labels <- character(8)
for (i in 1:8) {
  qvalue_labels[i] <- sprintf("[%d, %d]", qvalue_breaks[i], qvalue_breaks[i + 1])
}

qvalue_labels[8] <- sprintf(">=%d", qvalue_breaks[8])
# 在右下角添加图例
legend(
  "bottomright",                     
  legend = qvalue_labels,            
  fill = qvalue_colors,              
  title = "-log10(qvalue)",         
  border = "black",                   
  bty = "n",                        
  cex = 1.5,                       
)
################################################################
# Second track: Number of background genes
#circos.clear()

qvalue_range <- range(go_df$log10qvalue)
qvalue_breaks <- seq(qvalue_range[1], qvalue_range[2], length.out = 9)

qvalue_colors <- colorRampPalette(c("#F0C284", "#E3625D"))(8)
sector2 = go_df$number_of_genes
s2 = factor(sector2)
current_genes = go_df$log10qvalue[CELL_META$sector.numeric.index]

circos.track(s1, ylim = c(0, 1),
             bg.col = "white",
             bg.border = NA,
             panel.fun = function(x, y) {
               current_value = go_df$number_of_genes[CELL_META$sector.numeric.index]
               current_ratio = current_value / max(go_df$number_of_genes)
               current_qvalue = go_df$log10qvalue[CELL_META$sector.numeric.index]
               
               bin <- findInterval(current_qvalue, qvalue_breaks, rightmost.closed = TRUE)
               bin <- min(bin, 8)  
               current_color <- qvalue_colors[bin]
               circos.rect(0, 0, current_value, 1, 
                           col = current_color, 
                           border = NA)
               circos.text(CELL_META$xcenter, CELL_META$ycenter, 
                           labels = current_value, 
                           cex = 1.8,
                           col = "black",
                           facing = "inside",
                           niceFacing = TRUE
               )
             })
#############################################################
# Third track: Number of selected genes
circos.track(s1, ylim = c(0, 1),
             bg.col = "white",
             bg.border = NA,
             panel.fun = function(x,y){
               current_value1 = go_df$number_of_select[CELL_META$sector.numeric.index]
               circos.rect(0, 0, current_value1, 1, 
                           col = "#73A5A2", 
                           border = NA)
               circos.text(CELL_META$xcenter, 0.5, 
                           #current_genes,
                           labels = go_df$number_of_select[CELL_META$sector.numeric.index], 
                           cex = 1.8,  
                           col = "black", 
                           facing = "inside",  
                           niceFacing = TRUE)  
             })
#####################################################################
# Fourth to eighth tracks: RichFactor thresholds
#1
circos.par(track.height = 0.1)
richfactor_breaks <- seq(0, 0.1, length.out = 3)
current_rf <- go_df$richfactor[CELL_META$sector.numeric.index]
current_category <- go_df$Category[CELL_META$sector.numeric.index]
track_colors <- ifelse(go_df$richfactor > 0.5, 
                       category_colors[go_df$Category], 
                       "grey81")

circos.track(s1, ylim = c(0, 1),
             bg.col = track_colors,
             bg.border = "white",
             track.height = 0.1)

#3
circos.par(track.height = 0.1)
richfactor_breaks <- seq(0, 0.1, length.out = 3)
current_rf <- go_df$richfactor[CELL_META$sector.numeric.index]
current_category <- go_df$Category[CELL_META$sector.numeric.index]
track_colors <- ifelse(go_df$richfactor > 0.4, 
                       category_colors[go_df$Category], 
                       "grey81")

circos.track(s1, ylim = c(0, 1),
             bg.col = track_colors,
             bg.border = "white",
             track.height = 0.1)

#5
circos.par(track.height = 0.1)
richfactor_breaks <- seq(0, 0.1, length.out = 3)
current_rf <- go_df$richfactor[CELL_META$sector.numeric.index]
current_category <- go_df$Category[CELL_META$sector.numeric.index]
track_colors <- ifelse(go_df$richfactor > 0.3, 
                       category_colors[go_df$Category], 
                       "grey81")

circos.track(s1, ylim = c(0, 1),
             bg.col = track_colors,
             bg.border = "white",
             track.height = 0.1)

#7
circos.par(track.height = 0.1)
richfactor_breaks <- seq(0, 0.1, length.out = 3)
current_rf <- go_df$richfactor[CELL_META$sector.numeric.index]
current_category <- go_df$Category[CELL_META$sector.numeric.index]
track_colors <- ifelse(go_df$richfactor > 0.2, 
                       category_colors[go_df$Category], 
                       "grey81")

circos.track(s1, ylim = c(0, 1),
             bg.col = track_colors,
             bg.border = "white",
             track.height = 0.1)
#9
circos.par(track.height = 0.1)
richfactor_breaks <- seq(0, 0.1, length.out = 3)
current_rf <- go_df$richfactor[CELL_META$sector.numeric.index]
current_category <- go_df$Category[CELL_META$sector.numeric.index]
track_colors <- ifelse(go_df$richfactor > 0.1, 
                       category_colors[go_df$Category], 
                       "grey81")

circos.track(s1, ylim = c(0, 1),
             bg.col = track_colors,
             bg.border = "white",
             track.height = 0.1)

########################################################################
# Center legend
inside_color <- c(
  "Number of Genes" = "#E3625D",  
  "Number of Select" = "#73A5A2", 
  "RichFactor" = "#8074C8"  
)
circos.info()
#画图
legend(
  "center",                   
  legend = c("Number of Genes", "Number of Select", "RichFactor"),
  pch = c(15, 15, 17), 
  col = c("#E3625D", "#73A5A2", "#8074C8"),          
  title = NA,               
  border = c("black", "black", "black"),                      
  bty = "n",                         
  pt.cex = 1.6,
  cex = 1.7                          
)
####################End plotting
dev.off()
circos.clear()
####################################################################KEGG + KEGG Circos
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)


# Set working directory
path="C:/Users/29077/Desktop/20251203-enrichment"
setwd(path)

# Load gene list
genes_raw <- read.table("gene_intersection_482.txt", header = FALSE)$V1
gene_list <- genes_raw[genes_raw != "gene_intersection_482"]
length(gene_list)

# Convert Symbol to Entrez ID
gene_entrez <- mapIds(org.Hs.eg.db,
                      keys = gene_list,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")
gene_entrez <- na.omit(gene_entrez)  

# KEGG enrichment analysis
options(timeout = 360)
ekegg <- enrichKEGG(gene         = gene_entrez,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.xlsx(as.data.frame(ekegg), "KEGG_enrichment.xlsx")
###################################################################Prepare KEGG data
all_kegg <- read.xlsx("KEGG_enrichment.xlsx")
kegg_sorted <- all_kegg[order(all_kegg$qvalue),]

# Select top terms by category
human_diseases <- kegg_sorted[kegg_sorted$category == "Human Diseases", ]
human_diseases_top <- human_diseases[1:8, ]

environmental <- kegg_sorted[kegg_sorted$category == "Environmental Information Processing", ]
environmental_top <- environmental[1:6, ]

organismal <- kegg_sorted[kegg_sorted$category == "Organismal Systems", ]
organismal_top <- organismal[1:6, ]

selected_kegg <- rbind(human_diseases_top, environmental_top, organismal_top)

# Extract variables
kegg_top_ids <- selected_kegg$ID[1:20]
print(kegg_top_ids)

kegg_top_Category <- selected_kegg$category[1:20]
print(kegg_top_Category)

kegg_top_Subcategory <- selected_kegg$subcategory[1:20]
print(kegg_top_Subcategory)

kegg_top_qvalue <- selected_kegg$qvalue[1:20]
print(kegg_top_qvalue)
kegg_top_log10qvalue <- -log10(kegg_top_qvalue)
print(kegg_top_log10qvalue)

kegg_richfactor <- selected_kegg$RichFactor[1:20]
print(kegg_richfactor)

kegg_genes_number <- sapply(strsplit(selected_kegg$BgRatio[1:20], "/"), function(x) as.numeric(x[1]))
print(kegg_genes_number)

kegg_genes_select <- sapply(strsplit(selected_kegg$GeneRatio[1:20], "/"), function(x) as.numeric(x[1]))
print(kegg_genes_select)

kegg_description <- selected_kegg$Description[1:20]
print(kegg_description)

# Build data frame
kegg_df <- data.frame(
  Number = 1:length(kegg_top_ids),
  KEGG_ID = kegg_top_ids,
  log10qvalue = kegg_top_log10qvalue,
  number_of_genes = kegg_genes_number, 
  number_of_select = kegg_genes_select, 
  richfactor = kegg_richfactor,
  Category = rep(kegg_top_Category),
  subcategory = kegg_top_Subcategory,
  Description = kegg_description,
  stringsAsFactors = FALSE
)

print(kegg_df)
###############################################################################Special definition
kegg_df$Highlight <- ifelse(kegg_df$subcategory %in% c("Neurodegenerative disease", "Nervous system"), 
                            "#FF4C4C", "none")

#######################################################################
# Plot KEGG Circos
par(family = "serif", font = 1)  


tiff("KEGG_Enrichment_Circos.tiff", width = 23000, height = 23000, res = 1600, compression = "lzw")

windowsFonts(Times = windowsFont("Times New Roman"))
par(family = "Times")

#################################################
# Colors and scales
category_colors <- c("Human Diseases" = "#8BA8F5", "Organismal Systems" = "#A6A4CC", "Environmental Information Processing" = "#B88AF5", "NA" = "#777777", "Cellular Processes" = "#F1B66D")
highlight_color <- "#FF4C4C"
scale_range <- c(0, 500) 

# Reorder to highlight specific pathways
kegg_df <- kegg_df[order(kegg_df$Category), ]
highlight_idx <- which(kegg_df$Highlight == "#FF4C4C")
non_highlight_idx <- which(kegg_df$Highlight != "#FF4C4C")
kegg_df <- kegg_df[c(highlight_idx, non_highlight_idx), ]

circos.clear()
circos.par(
  start.degree = 270, 
  clock.wise = TRUE,
  track.margin = c(0, 0),
  "track.height" = 0.1
)
sector = kegg_df$KEGG_ID
s1 = factor(sector, levels = kegg_df$KEGG_ID)
circos.initialize(
  s1,
  xlim = scale_range
)
kegg_df$Category[is.na(kegg_df$Category)] <- "NA"
##############################################################
track_colors <- ifelse(kegg_df$Highlight == "#FF4C4C", 
                       highlight_color, 
                       category_colors[kegg_df$Category])
##########################################################################
# Outer track: KEGG ID / Description
circos.track(s1, ylim = c(0,1),
             bg.col = track_colors,
             panel.fun = function(x,y){
               current_desc <- kegg_df$Description[CELL_META$sector.numeric.index]
               display_text <- ifelse(current_desc == "Alzheimer disease", 
                                      "Alzheimer disease", 
                                      kegg_df$KEGG_ID[CELL_META$sector.numeric.index])
               text_col <- ifelse(kegg_df$Highlight[CELL_META$sector.numeric.index] == "#FF4C4C",
                                  "white", "black")
               circos.text(CELL_META$xcenter,
                           CELL_META$cell.ylim[1] + 0.7,
                           labels = display_text,  
                           facing = "bending.inside",
                           niceFacing = TRUE,
                           col = text_col,
                           cex = 2.0,
                           font = 1  
               )
               circos.axis(labels.cex = 1.3,
                           major.tick.length = 0.35,
                           major.at = c(0, 100, 300, 500),
                           labels.font = 1
               )
             })
##################################################################
# Second track: Number of background genes
qvalue_range <- range(kegg_df$log10qvalue)
qvalue_breaks <- seq(qvalue_range[1], qvalue_range[2], length.out = 9)

qvalue_colors <- colorRampPalette(c("navajowhite", "firebrick1"))(8)
sector2 = kegg_df$log10qvalue
s2 = factor(sector2)
current_genes = kegg_df$log10qvalue[CELL_META$sector.numeric.index]
circos.track(s1, ylim = c(0, 1),
             bg.col = "white",
             bg.border = NA,
             panel.fun = function(x, y) {
               current_value = kegg_df$number_of_genes[CELL_META$sector.numeric.index]
               current_ratio = current_value / max(kegg_df$number_of_genes)
               current_qvalue = kegg_df$log10qvalue[CELL_META$sector.numeric.index]
               
               bin <- findInterval(current_qvalue, qvalue_breaks, rightmost.closed = TRUE)
               bin <- min(bin, 8)  
               current_color <- qvalue_colors[bin]
               circos.rect(0, 0, current_value, 1, 
                           col = current_color, 
                           border = NA)
               circos.text(CELL_META$xcenter, CELL_META$ycenter, 
                           labels = current_value, 
                           cex = 1.8,
                           col = "black",
                           facing = "inside",
                           niceFacing = TRUE,
                           font = 1
               )
             })
#####################################################################
# Third track: Number of selected genes
circos.track(s1, ylim = c(0, 1),
             bg.col = "white",
             bg.border = NA,
             panel.fun = function(x,y){
               current_value1 = kegg_df$number_of_select[CELL_META$sector.numeric.index]
               circos.rect(0, 0, current_value1, 1, 
                           col = "#73A5A2", 
                           border = NA)
               circos.text(CELL_META$xcenter, 0.5, 
                           labels = kegg_df$number_of_select[CELL_META$sector.numeric.index], 
                           cex = 1.8,
                           col = "black",
                           facing = "inside",
                           niceFacing = TRUE,
                           font = 1)
             })
#####################################################################
# Inner tracks: RichFactor thresholds
track_colors_threshold <- function(threshold) {
  ifelse(
    kegg_df$richfactor > threshold,
    ifelse(
      kegg_df$Highlight == "#FF4C4C",
      highlight_color,            
      category_colors[kegg_df$Category]  
    ),
    "grey81" 
  )
}

# 1
circos.par(track.height = 0.1)
circos.track(
  s1, ylim = c(0, 1),
  bg.col = track_colors_threshold(0.5),
  bg.border = "white",
  track.height = 0.1
)

#3
circos.par(track.height = 0.1)
circos.track(
  s1, ylim = c(0, 1),
  bg.col = track_colors_threshold(0.4),
  bg.border = "white",
  track.height = 0.1
)

#5
circos.par(track.height = 0.1)
circos.track(
  s1, ylim = c(0, 1),
  bg.col = track_colors_threshold(0.3),
  bg.border = "white",
  track.height = 0.1
)

#7
circos.par(track.height = 0.1)
circos.track(
  s1, ylim = c(0, 1),
  bg.col = track_colors_threshold(0.2),
  bg.border = "white",
  track.height = 0.1
)

#9
circos.par(track.height = 0.1)
circos.track(
  s1, ylim = c(0, 1),
  bg.col = track_colors_threshold(0.1),
  bg.border = "white",
  track.height = 0.1
)

#########################################################
# Legend: KEGG categories (bottom-left)
category_color <- c(
  "Human Diseases" = "#8BA8F5",  
  "Organismal Systems" = "#A6A4CC", 
  "Environmental Information Processing" = "#B88AF5"
)

legend(
  "bottomleft",                      
  legend = names(category_color),   
  fill = category_color,           
  title = "CATEGORY",                
  border = NA,                   
  bty = "n",                        
  cex = 1.4,                        
  text.font = 1,                   
  title.font = 2                 
)
#################################################################
# Legend: -log10(qvalue) gradient (bottom-right)
qvalue_range <- range(kegg_df$log10qvalue)
qvalue_breaks <- round(seq(qvalue_range[1], qvalue_range[2], length.out = 9))
qvalue_colors <- colorRampPalette(c("navajowhite", "firebrick1"))(8) 
qvalue_labels <- character(8)
for (i in 1:8) {
  qvalue_labels[i] <- sprintf("[%d, %d]", qvalue_breaks[i], qvalue_breaks[i + 1])
}

qvalue_labels[8] <- sprintf(">=%d", qvalue_breaks[8])

legend(
  "bottomright",                     
  legend = qvalue_labels,          
  fill = qvalue_colors,             
  title = "-log10(qvalue)",         
  border = "black",                
  bty = "n",                       
  cex = 1.5,                         
  text.font = 1,                   
  title.font = 2                  
)
#################################################
# Center legend
legend(
  "center",                      
  legend = c("Number of Genes", "Number of Select", "RichFactor"),
  pch = c(15, 15, 17), 
  col = c("navajowhite", "#73A5A2", "#456ade"),          
  title = NA,                
  border = c("black", "black", "black"),                      
  bty = "n",                         
  pt.cex = 1.6,
  cex = 1.7,                          
  text.font = 1                    
)
#####################################################
dev.off()

#########################################################################
# Dot plots
#######################################################################GO dot plot

library(readxl)
library(dplyr)
library(ggplot2)
library(extrafont)


go_data <- read.xlsx("top10_GO.xlsx", sheet = 1)


go_data$GeneRatio <- as.numeric(sapply(strsplit(go_data$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
go_data$RichFactor <- as.numeric(go_data$RichFactor)
go_data$qvalue <- as.numeric(go_data$qvalue)
go_data$Count <- as.numeric(go_data$Count)


p_go <- ggplot(go_data, aes(x = RichFactor, y = reorder(Description, RichFactor))) +
  geom_point(aes(size = Count, fill = -log10(qvalue)), 
             shape = 21,
             color = "white",
             stroke = 0) +
  scale_fill_gradient(low = "#A8CBDF", high = "#7895C1", 
                      name = "-log10(qvalue)") +
  scale_size_continuous(name = "Gene Count", 
                        range = c(4, 12),
                        guide = guide_legend(override.aes = list(fill = "black", color = "black"))) +
  labs(x = "Rich Factor", 
       y = "GO Term",
       title = "GO Enrichment") +
  theme_bw() +
  theme(

    text = element_text(family = "Times New Roman", size = 12),
    

    axis.text.y = element_text(size = 18, color = "black", family = "Times New Roman"),
    axis.text.x = element_text(size = 18, color = "black", family = "Times New Roman"),
    

    axis.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title.x = element_text(margin = margin(t = 10), family = "Times New Roman"),
    axis.title.y = element_text(margin = margin(r = 10), family = "Times New Roman"),
    

    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                              margin = margin(b = 15), family = "Times New Roman"),
    

    legend.title = element_text(family = "Times New Roman", size = 13, face = "bold"),
    legend.text = element_text(family = "Times New Roman", size = 12)
  )


print(p_go)

#save
ggsave("GO-top10-dot.tiff", 
       plot = p_go,
       device = "tiff",
       width = 10, 
       height = 10,
       dpi = 800,
       compression = "lzw")
####################################################KEGG

library(openxlsx)
library(ggplot2)
library(stringr)
library(extrafont)



data <- read.xlsx("top10_KEGG.xlsx", sheet = 1)

data$GeneRatio <- as.numeric(sapply(strsplit(data$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
data$RichFactor <- as.numeric(data$RichFactor)
data$qvalue <- as.numeric(data$qvalue)
data$Count <- as.numeric(data$Count)


p_black_legend <- ggplot(data, aes(x = RichFactor, y = reorder(Description, RichFactor))) +
  geom_point(aes(size = Count, fill = -log10(qvalue)), 
             shape = 21,
             color = "white",
             stroke = 0) +
  scale_fill_gradient(low = "#A6A4CC", high = "#B88AF5", 
                      name = "-log10(qvalue)") +
  scale_size_continuous(name = "Gene Count", 
                        range = c(4, 12),
                        guide = guide_legend(override.aes = list(fill = "black", color = "black"))) +
  labs(x = "Rich Factor", 
       y = "KEGG Pathway",
       title = "KEGG Pathways Enrichment") +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    

    axis.text.y = element_text(size = 18, color = "black", family = "Times New Roman"),
    axis.text.x = element_text(size = 18, color = "black", family = "Times New Roman"),
    
    axis.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title.x = element_text(margin = margin(t = 10), family = "Times New Roman"),
    axis.title.y = element_text(margin = margin(r = 10), family = "Times New Roman"),

    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                              margin = margin(b = 15), family = "Times New Roman"),
    
    legend.title = element_text(family = "Times New Roman", size = 13, face = "bold"),
    legend.text = element_text(family = "Times New Roman", size = 12)
  )


print(p_black_legend)


ggsave("KEGG-top10-dot.tiff", 
       plot = p_black_legend,
       device = "tiff",
       width = 10, 
       height = 10,
       dpi = 800,
       compression = "lzw")

