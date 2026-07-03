setwd("/home/xiaoning/Documents/projects/MD-sun/PFOS-CKM/data/GSE43292")

library(GEOquery)

gset <- getGEO("GSE43292", GSEMatrix = TRUE, AnnotGPL = FALSE)[[1]]


expr_mat <- exprs(gset)


write.table(
  expr_mat,
  file = "GSE43292_expression_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

dim(expr_mat)  


expr_mat <- read.table(
  "GSE43292_expression_matrix.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)


coldata <- read.table(
  "GSE43292_deg_sampleinfo.txt",
  header = TRUE,
  sep = "\t"
)


coldata <- coldata[coldata$SampleID %in% colnames(expr_mat), ]
coldata <- coldata[match(colnames(expr_mat), coldata$SampleID), ]

stopifnot(nrow(coldata) == 64)
stopifnot(all(colnames(expr_mat) == coldata$SampleID))


gpl <- getGEO("GPL6244")
anno <- Table(gpl)

get_symbol <- function(x){
  if (is.na(x) || x == "---") return(NA)
  first <- strsplit(x, " /// ")[[1]][1]
  parts <- strsplit(first, " // ")[[1]]
  if (length(parts) >= 2) return(trimws(parts[2]))
  return(NA)
}

anno$Gene <- sapply(anno$gene_assignment, get_symbol)
anno2 <- anno[, c("ID", "Gene")]


expr_mat$ID <- rownames(expr_mat)
merged <- merge(anno2, expr_mat, by = "ID")

merged <- merged[!is.na(merged$Gene), ]
merged <- merged[merged$Gene != "" & merged$Gene != "---", ]

library(dplyr)

expr_gene <- merged %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), mean))


write.csv(expr_gene, "expression.csv", row.names = FALSE)


coldata$tissue <- factor(
  coldata$tissue,
  levels = c("Macroscopically intact tissue", "Atheroma plaque")
)


expr_mat$ID <- NULL

label <- data.frame(
  Sample = colnames(expr_mat),
  Group = ifelse(coldata$tissue == "Atheroma plaque", 1, 0)
)

write.csv(label, "label.csv", row.names = FALSE)

library(dplyr)


expr <- read.csv("expression.csv", row.names = 1, check.names = FALSE)
label <- read.csv("label.csv")


expr$ID <- NULL

active_genes <- c("DES", "SPP1", "ELN", "FABP4", "MMP7")


t_test_results <- sapply(active_genes, function(gene) {
  if (!gene %in% rownames(expr)) return(NULL)
  
  values <- as.numeric(expr[gene, ])
  
  group0 <- values[label$Group == 0]
  group1 <- values[label$Group == 1]
  
  if (length(group0) < 2 || length(group1) < 2) return(NULL)
  
  test <- t.test(group1, group0)
  
  data.frame(
    Gene = gene,
    mean_plaque  = mean(group1, na.rm = TRUE),
    mean_intact = mean(group0, na.rm = TRUE),
    log2FC       = mean(group1, na.rm = TRUE) - mean(group0, na.rm = TRUE),
    p_value      = test$p.value
  )
}, simplify = FALSE)

t_test_results <- bind_rows(t_test_results)


t_test_results <- t_test_results %>%
  mutate(
    signif_star = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE              ~ "ns"
    )
  )

print(t_test_results)



library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggthemes)



expr <- read.csv("expression.csv", row.names = 1, check.names = FALSE)
expr$ID <- NULL

label <- read.csv("label.csv")


active_genes <- c("DES", "SPP1", "ELN", "FABP4", "MMP7")


t_test_results <- sapply(active_genes, function(gene) {
  if (!gene %in% rownames(expr)) return(NULL)
  
  values <- as.numeric(expr[gene, ])
  
  group0 <- values[label$Group == 0]
  group1 <- values[label$Group == 1]
  
  if (length(group0) < 2 || length(group1) < 2) return(NULL)
  
  test <- t.test(group1, group0)
  
  data.frame(
    Gene = gene,
    mean_plaque  = mean(group1),
    mean_intact = mean(group0),
    log2FC       = mean(group1) - mean(group0),
    p_value      = test$p.value
  )
}, simplify = FALSE)

t_test_results <- bind_rows(t_test_results) %>%
  mutate(
    signif_star = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE              ~ "ns"
    )
  )

print(t_test_results)


expr5 <- expr[rownames(expr) %in% active_genes, ]

expr_violin_long <- expr5 %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "Sample", values_to = "expression") %>%
  left_join(label, by = c("Sample" = "Sample"))


y_positions <- expr_violin_long %>%
  group_by(gene) %>%
  summarise(max_expr = max(expression, na.rm = TRUE)) %>%
  mutate(y_position = max_expr + 0.5)

signif_data <- left_join(y_positions, t_test_results, by = c("gene" = "Gene")) %>%
  filter(signif_star != "ns")


p <- ggplot(expr_violin_long,
            aes(x = gene, y = expression, fill = factor(Group))) +
  geom_violin(trim = FALSE, color = "white",
              position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9),
               outlier.shape = NA) +
  scale_fill_manual(
    values = c("0" = "lightskyblue1", "1" = "lightgoldenrod1"),
    labels = c("0" = "Intact", "1" = "Plaque"),
    name  = "Tissue"
  ) +
  

  labs(
    title = "Validation of Candidate Genes in GSE43292",
    x = "Gene",
    y = "log2(Expression)"
  ) +
  
  theme_few(base_family = "Times New Roman") +
  
  theme(

    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, colour = "black"),
    

    axis.title.x = element_text(size = 14, colour = "black"),
    axis.text.x  = element_text(size = 13, angle = 0, hjust = 0.5, colour = "black"),
    

    axis.title.y = element_text(size = 14, colour = "black"),
    axis.text.y  = element_text(size = 12, colour = "black"),
    

    legend.position      = c(0.05, 0.05),
    legend.justification = c(0, 0),
    legend.title  = element_text(size = 14, face = "bold", colour = "black"),
    legend.text   = element_text(size = 13),
    legend.key    = element_blank(),
    legend.background = element_blank()
  )


if (nrow(signif_data) > 0) {
  for (i in 1:nrow(signif_data)) {
    gene_name <- signif_data$gene[i]
    x_pos <- which(levels(factor(expr_violin_long$gene)) == gene_name)
    y_pos <- signif_data$max_expr[i] + 0.8
    
    x_intact <- x_pos - 0.2
    x_plaque <- x_pos + 0.2
    
    p <- p +
      geom_segment(x = x_intact,  xend = x_intact,
                   y = y_pos - 0.1, yend = y_pos,
                   color = "black", linewidth = 0.8) +
      geom_segment(x = x_plaque, xend = x_plaque,
                   y = y_pos - 0.1, yend = y_pos,
                   color = "black", linewidth = 0.8) +
      geom_segment(x = x_intact, xend = x_plaque,
                   y = y_pos, yend = y_pos,
                   color = "black", linewidth = 0.8) +
      annotate("text",
               x = (x_intact + x_plaque) / 2,
               y = y_pos + 0.3,
               label = signif_data$signif_star[i],
               size = 6, color = "black", fontface = "bold")
  }
}

print(p)


num_genes <- length(unique(expr_violin_long$gene))
plot_width <- max(10, num_genes * 0.8)  


# ggsave("violin_plot_genes_20260630.pdf", plot = p, width = plot_width, height = 6,
#        family = "Times")

ggsave(
  "violin_plot_genes_20260630.png",
  plot  = p,
  width  = plot_width,
  height = 6,
  dpi    = 300
)


