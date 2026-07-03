
setwd("/home/xiaoning/Documents/projects/MD-sun/PFOS-CKM/data/GSE120521")


##########################################################################################################################

BiocManager::install("glmnet")

BiocManager::install("readr")
BiocManager::install("readxl")
BiocManager::install("dplyr")
BiocManager::install("sva")
BiocManager::install("edgeR")
# ------------------------------

library(DESeq2)    
library(glmnet)    
library(edgeR)     
library(sva)      


library(DESeq2)    
library(glmnet)    
library(readxl)    
library(tidyverse) 
library(openxlsx) 
library(sva)       

library(org.Hs.eg.db)
library(AnnotationDbi)



expr <- read.table(
  "GSE120521_raw_counts_GRCh38.p13_NCBI.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)


genes <- read.table(
  "intersection_genes_stabilityofregion_483.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)[,1]

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys       = genes,
  keytype    = "SYMBOL",
  columns    = c("ENTREZID")
)
gene_map <- gene_map[!is.na(gene_map$ENTREZID), ]
entrez_only <- gene_map$ENTREZID
write.table(
  entrez_only,
  "intersection_genes_stabilityofregion_483_EntrezID.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
deg <- read.table("intersection_genes_stabilityofregion_483_EntrezID.txt", header = FALSE)
colnames(deg) <- "gene"


sampleLabel <- read.table(
  "GSE120521_deg_sampleinfo.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1
)

head(sampleLabel)


expr_sub <- expr[rownames(expr) %in% deg$gene, ]


cat("gene", nrow(expr_sub), "\n")

colData <- data.frame(row.names = colnames(expr_sub),
                      condition = factor(sampleLabel$stabilityofregion))


dds <- DESeqDataSetFromMatrix(countData = round(expr_sub),
                              colData = colData,
                              design = ~ condition)


dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)


expr_log <- log2(norm_counts + 1)


expr_scaled <- t(scale(t(expr_log))) 


X <- t(expr_scaled)   
y <- factor(sampleLabel$stabilityofregion)  


stopifnot(nrow(X) == length(y))


set.seed(123)  
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5)


cat("best lambda.min:", cvfit$lambda.min, "\n")
cat("best lambda.1se:", cvfit$lambda.1se, "\n")


coef_min <- coef(cvfit, s = "lambda.min")
active_index <- which(coef_min != 0)
active_genes <- rownames(coef_min)[active_index]
active_genes <- active_genes[active_genes != "(Intercept)"]

cat("catch gene：\n")
print(active_genes)

write.csv(active_genes, file = "LASSO_selected_genes_symbol.csv", row.names = FALSE)
mapped_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = active_genes,      
  columns = "SYMBOL",     
  keytype = "ENTREZID"   
)
head(mapped_genes)
genes_only <- mapped_genes$SYMBOL
write.table(
  genes_only,
  "LASSO_selected_genes.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)



library(glmnet)

pdf("LASSO_CV_and_CoefficientPaths.pdf", width = 10, height = 5,family = "Times")
par(mfrow = c(1,2))


plot(cvfit,
     xlab = "Log(lambda)",
     ylab = "Mean-Squared Error")


plot(cvfit$glmnet.fit,
     xvar = "lambda",
     label = TRUE,
     xlab = "Log(lambda)",
     ylab = "Coefficients")
dev.off()

cat("Save PDF：LASSO_CV_and_CoefficientPaths.pdf\n")



coef_df <- as.data.frame(as.matrix(coef_min))
colnames(coef_df) <- "Coefficient"
coef_df <- coef_df[coef_df$Coefficient != 0, , drop = FALSE]
write.csv(coef_df, "LASSO_selected_genes.csv", quote = FALSE)

cat("Save LASSO_selected_genes.csv\n")


coef(cvfit, s = "lambda.min")

library(ggplot2)

expr_log <- log2(expr + 1)
expr7 <- expr_log[rownames(expr_log) %in% c("6696","2006","2167"), ]

library(dplyr)
library(tidyr)
library(tibble)

sampleLabel2 <- sampleLabel %>%
  rownames_to_column(var = "sample")  

expr7_long <- as.data.frame(t(expr7)) %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(-sample, names_to = "gene", values_to = "expression") %>%
  left_join(sampleLabel2, by = "sample")  

ggplot(expr7_long, aes(x = factor(stabilityofregion), y = expression, fill = factor(stabilityofregion))) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") +
  theme_bw() +
  labs(x = "Group", y = "log2(Normalized Counts + 1)")




#install.packages("dplyr")
#install.packages("rstatix")
#install.packages("ggsignif")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggthemes)
library(rstatix)
library(ggsignif)
library(dplyr)
# -----------------------------

active_genes <- c("6696","2006","2167","1674","4316") #"12"

mapped_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = active_genes,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)


entrez2symbol <- setNames(mapped_genes$SYMBOL, mapped_genes$ENTREZID)

# -----------------------------

expr_log <- log2(expr + 1)
expr_violin <- expr_log[rownames(expr_log) %in% active_genes, ]


rownames(expr_violin) <- entrez2symbol[rownames(expr_violin)]

# -----------------------------

sampleLabel2 <- sampleLabel %>%
  rownames_to_column(var = "sample") 

expr_violin_long <- as.data.frame(t(expr_violin)) %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(-sample, names_to = "gene", values_to = "expression") %>%
  left_join(sampleLabel2, by = "sample")  


deg_symbol <- read.table("PTC_vs_ATC_DEG_with_symbol.txt", 
                  sep = "\t", header = TRUE, check.names = FALSE)

deg2 <- deg_symbol[, c("SYMBOL", "padj")]

deg2 <- deg2 %>% filter(SYMBOL %in% rownames(expr_violin))

deg2 <- deg2 %>%
  mutate(signif_star = case_when(
    padj <= 0.001 ~ "***",
    padj <= 0.01  ~ "**",
    padj <= 0.05  ~ "*",
    TRUE ~ "ns"
  ))

y_positions <- expr_violin_long %>%
  group_by(gene) %>%
  summarise(
    max_expr = max(expression, na.rm = TRUE),
    min_expr = min(expression, na.rm = TRUE)
  ) %>%
  mutate(
    y_position = max_expr + 0.5 
  )

signif_data <- left_join(y_positions, deg2, 
                         by = c("gene" = "SYMBOL")) %>%
  filter(signif_star != "ns")      


expr_violin_long$stabilityofregion <- factor(
  expr_violin_long$stabilityofregion,
  levels = c("unstable", "stable")
)
p <- ggplot(expr_violin_long, aes(x = gene, y = expression, fill = stabilityofregion)) +
  geom_violin(trim = FALSE, color = "white", position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("stable" = "lightskyblue1", "unstable" = "lightgoldenrod1")) +
  theme_few() +
  labs(x = "Gene", y = "log2(Normalized Counts + 1)", fill = "Stability Region") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    
    legend.position = c(0.05, 0.05),   
    legend.justification = c(0, 0),    
    
    legend.background = element_blank(),  
    legend.key = element_blank()          
  )


if (nrow(signif_data) > 0) {
  for (i in 1:nrow(signif_data)) {
    gene_name <- signif_data$gene[i]
    x_pos <- which(levels(factor(expr_violin_long$gene)) == gene_name)
    y_pos <- signif_data$max_expr[i] + 1.0 

    x_stable <- x_pos - 0.2
    x_unstable <- x_pos + 0.2

    p <- p +

      geom_segment(x = x_stable, xend = x_stable,
                   y = y_pos - 0.1, yend = y_pos,
                   color = "black", linewidth = 0.8) +

      geom_segment(x = x_unstable, xend = x_unstable,
                   y = y_pos - 0.1, yend = y_pos,
                   color = "black", linewidth = 0.8) +

      geom_segment(x = x_stable, xend = x_unstable,
                   y = y_pos, yend = y_pos,
                   color = "black", linewidth = 0.8) +

      annotate("text",
               x = (x_stable + x_unstable) / 2, 
               y = y_pos + 0.3,
               label = signif_data$signif_star[i],   
               size = 6, color = "black", fontface = "bold")
  }
}

print(p)


num_genes <- length(unique(expr_violin_long$gene))
plot_width <- max(10, num_genes * 0.8)  


ggsave("violin_plot_genes_20260630.pdf", plot = p, width = plot_width, height = 6,
       family = "Times")


###########################################################################################################################
##############################################################
BiocManager::install("caret")
BiocManager::install("e1071")



library(e1071)      
library(caret)      
library(tidyverse)  
library(readxl)     


expr <- read.table(
  "GSE120521_raw_counts_GRCh38.p13_NCBI.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)


genes <- read.table(
  "intersection_genes_stabilityofregion_483.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)[,1]

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys       = genes,
  keytype    = "SYMBOL",
  columns    = c("ENTREZID")
)
gene_map <- gene_map[!is.na(gene_map$ENTREZID), ]
entrez_only <- gene_map$ENTREZID
write.table(
  entrez_only,
  "intersection_genes_stabilityofregion_483_EntrezID.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
deg <- read.table("intersection_genes_stabilityofregion_483_EntrezID.txt", header = FALSE)
colnames(deg) <- "gene"




sampleLabel <- read.table(
  "GSE120521_deg_sampleinfo.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1
)

head(sampleLabel)



expr_sub <- expr[rownames(expr) %in% deg$gene, ]


cat("catch gene:", nrow(expr_sub), "\n")


colData <- data.frame(row.names = colnames(expr_sub),
                      condition = factor(sampleLabel$stabilityofregion))

dds <- DESeqDataSetFromMatrix(countData = round(expr_sub),
                              colData = colData,
                              design = ~ condition)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

expr_log <- log2(norm_counts + 1)



expr17 <- expr_log[rownames(expr_log) %in% deg$gene, ]

cat("extract gene:", nrow(expr17), "\n")


expr_scaled <- t(scale(t(expr17)))  


X <- t(expr_scaled)  



sampleLabel <- read.table(
  "GSE120521_deg_sampleinfo.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1
)



X <- X[match(sampleLabel$SampleID, rownames(X)), ]
y <- factor(sampleLabel$stabilityofregion)

stopifnot(rownames(X) == sampleLabel$SampleID)


set.seed(123)
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]


svm_model <- svm(X_train, y_train,
                 kernel = "linear",
                 scale = FALSE,
                 probability = TRUE)


pred <- predict(svm_model, X_test)
conf_matrix <- confusionMatrix(pred, y_test)

print(conf_matrix)


library(pROC)


pred_prob <- attr(predict(svm_model, X_test, probability = TRUE), "probabilities")[, 2]

roc_obj <- roc(as.numeric(y_test), pred_prob)
plot(roc_obj, print.auc = TRUE, main = "SVM ROC Curve")


saveRDS(svm_model, "svm_model.rds")
write.csv(as.data.frame(conf_matrix$table), "svm_confusion_matrix.csv")

cat("save SVM model svm_model.rds，save matrix svm_confusion_matrix.csv\n")


##############################################################


library(caret)
library(ggplot2)
library(cowplot)   
library(e1071)     


data(iris)
X <- iris[, 1:4]
y <- iris$Species


set.seed(123)

ctrl <- rfeControl(functions = caretFuncs,
                   method = "cv",
                   number = 5)


svmProfile <- rfe(x = X,
                  y = y,
                  sizes = 1:4,
                  rfeControl = ctrl,
                  method = "svmLinear")

print(svmProfile)


results <- svmProfile$results
results


results$CVError <- 1 - results$Accuracy

best_idx <- which.max(results$Accuracy)
best_n   <- results$Variables[best_idx]
best_acc <- round(results$Accuracy[best_idx], 3)
best_err <- round(results$CVError[best_idx], 3)


p1 <- ggplot(results, aes(x = Variables, y = Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_text(aes(x = best_n, y = best_acc, 
                label = paste0(best_n, " - ", best_acc)),
            color = "red", vjust = -1) +
  theme_bw() +
  labs(x = "Number of Features", y = "5 × CV Accuracy")


p2 <- ggplot(results, aes(x = Variables, y = CVError)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_text(aes(x = best_n, y = best_err, 
                label = paste0(best_n, " - ", best_err)),
            color = "red", vjust = -1) +
  theme_bw() +
  labs(x = "Number of Features", y = "5 × CV Error")


pdf("SVM_RFE_Curves.pdf", width = 10, height = 5)
plot_grid(p1, p2, ncol = 2, labels = c("A", "B"))
dev.off()



svm_weights <- coef(svm_model)  

abs_weights <- abs(svm_weights)
top_genes <- names(sort(abs_weights, decreasing = TRUE))[1:4]  
print(top_genes) 

write.table(top_genes, file = "SVMtopgenes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)




BiocManager::install("randomForest")
BiocManager::install("e1071")
BiocManager::install("readr")

library(DESeq2)        
library(randomForest)  
library(caret)        
library(ggplot2)       
library(readxl)       
library(tidyverse)    


expr <- read.table(
  "GSE120521_raw_counts_GRCh38.p13_NCBI.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)



sampleLabel <- read.table(
  "GSE120521_deg_sampleinfo.txt",
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1
)

head(sampleLabel)

common_samples <- intersect(colnames(expr), rownames(sampleLabel))


expr <- expr[, common_samples, drop = FALSE]
sampleLabel <- sampleLabel[common_samples, , drop = FALSE]



genes <- read.table(
  "intersection_genes_stabilityofregion_483.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)[,1]

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys       = genes,
  keytype    = "SYMBOL",
  columns    = c("ENTREZID")
)
gene_map <- gene_map[!is.na(gene_map$ENTREZID), ]
entrez_only <- gene_map$ENTREZID
write.table(
  entrez_only,
  "intersection_genes_stabilityofregion_483_EntrezID.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
deg <- read.table("intersection_genes_stabilityofregion_483_EntrezID.txt", header = FALSE)
colnames(deg) <- "gene"


common_genes <- intersect(rownames(expr), entrez_only)
expr_selected <- expr[common_genes, , drop = FALSE]


any(is.na(expr_selected))  
expr_selected[is.na(expr_selected)] <- apply(expr_selected, 2, function(x) mean(x, na.rm = TRUE))  


sampleLabel <- as.data.frame(sampleLabel)
expr_selected <- as.matrix(expr_selected)
mode(expr_selected)  

col_data <- data.frame(row.names = colnames(expr_selected), condition = factor(sampleLabel$stabilityofregion))
dds <- DESeqDataSetFromMatrix(countData = expr_selected, colData = col_data, design = ~ condition)


dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)
expr_log <- log2(norm_counts + 1)


y <- factor(sampleLabel$stabilityofregion)


accuracy_results <- data.frame(Num_Features = integer(0), Accuracy = numeric(0))


num_genes <- nrow(expr_log)

expr_log_t <- t(expr_log) 


dim(expr_log_t)



for (num_features in 1:num_genes) {

  X <- expr_log_t[, 1:num_features, drop = FALSE]
  

  rf_model <- randomForest(X, y, importance = TRUE, ntree = 500)
  

  accuracy <- 1 - tail(rf_model$err.rate[, 1], 1) 
  

  accuracy_results <- rbind(accuracy_results, data.frame(Num_Features = num_features, Accuracy = accuracy))
}


library(ggplot2)
p <- ggplot(accuracy_results, aes(x = Num_Features, y = Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  theme_minimal() +
  labs(x = "Number of Features", y = "Accuracy (Cross-Validation)") +
  theme(text = element_text(size = 14)) +
  scale_x_continuous(limits = c(0, 10))

ggsave("RandomForest_Accuracy_Curve.pdf", plot = p, width = 6, height = 5, family = "Times")

cat("Save PDF：RandomForest_Accuracy_Curve.pdf\n")


library(AnnotationDbi)
library(org.Hs.eg.db)


rf_model_final <- randomForest(expr_log_t, y, importance = TRUE, ntree = 500)


feature_importance <- importance(rf_model_final)[, 1]


top_indices <- order(feature_importance, decreasing = TRUE)[1:3]
top_values  <- feature_importance[top_indices]


top_gene_ids <- colnames(expr_log_t)[top_indices]

print(data.frame(
  GeneID = top_gene_ids,
  Importance = top_values
))

entrez_ids <- as.character(top_gene_ids)

mapped_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = entrez_ids,
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
)

mapped_genes


write.csv(accuracy_results, "random_forest_accuracy_results.csv", row.names = FALSE)

write.csv(mapped_genes, "random_foresttop_important_genes.csv", row.names = TRUE)




library(VennDiagram)
library(readr)
library(tidyverse)



rf_genes <- read.csv("random_foresttop_important_genes.csv", row.names = 1)
rf_genes <- as.character(rf_genes[,1])

lasso_genes <- read.csv("LASSO_selected_genes.csv", row.names = 1)
lasso_genes <- rownames(lasso_genes)

lasso_data <- read.csv("LASSO_selected_genes.csv", row.names = 1)

lasso_data <- lasso_data[!rownames(lasso_data) %in% "(Intercept)", , drop = FALSE]

lasso_genes <- rownames(lasso_data[lasso_data$Coefficient != 0, , drop = FALSE])


svm_genes <- read.table("SVMtopgenes.txt", header = FALSE)$V1


venn.plot <- venn.diagram(
  x = list(
    LASSO = lasso_genes,
    RandomForest = rf_genes
  ),
  filename = "Venn_RF_LASSO_SVM.tiff",
  fill = c("#E5C1B3", "#F1D27B"),
  alpha = 0.5,
  scaled = FALSE,
  euler.d = FALSE,
  rotation.degree = 0,
  lwd = c(0, 0),

  cex = 1.2,
  fontface = "bold",
  fontfamily = "Times New Roman",
  cat.cex = 1.1,
  cat.fontfamily = "Times New Roman",
  cat.fontface = "bold",
  cat.col = c("#C76A49", "#F1A700"),
  

  cat.pos = c(-6, 6),
  cat.dist = c(0.09, 0.09),
  cat.default.pos = "outer",
  margin = 0.20
) 

intersection_genes <- Reduce(intersect, list(rf_genes, lasso_genes, svm_genes))

print(intersection_genes)


write.table(intersection_genes, file = "Intersection_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)





install.packages(
  "graphlayouts",
  dependencies = TRUE,
  repos = "https://cloud.r-project.org"
)
library(graphlayouts)
install.packages("units", dependencies = TRUE)
library(units)
install.packages("sf", dependencies = TRUE)
library(sf)

install.packages(
  "ggraph",
  dependencies = TRUE,
  repos = "https://cloud.r-project.org"
)
library(ggraph)
library(tidyverse)
library(igraph)




nodes <- tibble(
  name = c(
    "Lasso+GBM", "Hubgene", "MCODE",
    "CCL2", "CXCL8", "PTGS2", "IL6", "TNFAIP3",
    "STAT1", "IRF7", "CXCL10", "ISG15"
  ),
  group = c(
    "Method","Method","Method",
    "Lasso+GBM","Lasso+GBM","Lasso+GBM","Lasso+GBM","Lasso+GBM",
    "MCODE","MCODE","MCODE","MCODE"
  )
)

edges <- tribble(
  ~from,        ~to,
  "Lasso+GBM",  "CCL2",
  "Lasso+GBM",  "CXCL8",
  "Lasso+GBM",  "PTGS2",
  "Lasso+GBM",  "IL6",
  "Lasso+GBM",  "TNFAIP3",
  "MCODE",      "STAT1",
  "MCODE",      "IRF7",
  "MCODE",      "CXCL10",
  "MCODE",      "ISG15",
  "Hubgene",    "IL6",
  "Hubgene",    "STAT1"
)

g <- graph_from_data_frame(
  d = edges,
  vertices = nodes,
  directed = FALSE
)

ggraph(g, layout = "circle") +
  geom_edge_arc(
    aes(color = after_stat(index)),
    curvature = 0.2,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(color = group, size = group),
    show.legend = FALSE
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3
  ) +
  scale_size_manual(
    values = c(
      "Method" = 8,
      "Lasso+GBM" = 4,
      "MCODE" = 4
    )
  ) +
  scale_color_manual(
    values = c(
      "Method" = "#4D4D4D",
      "Lasso+GBM" = "#6BAED6",
      "MCODE" = "#E6550D"
    )
  ) +
  theme_void()












