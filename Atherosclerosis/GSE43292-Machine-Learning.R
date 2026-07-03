
setwd("/home/xiaoning/Documents/projects/MD-sun/PFOS-CKM/data/GSE43292新")

library(dplyr)
library(caret)
library(glmnet)
library(e1071)
library(randomForest)
library(pROC)
library(ggplot2)


expr <- read.csv("expression.csv", row.names = 1, check.names = FALSE)
expr$ID <- NULL

label <- read.csv("label.csv")


active_genes <- c("DES", "SPP1", "ELN", "FABP4", "MMP7")
expr5 <- expr[rownames(expr) %in% active_genes, ]


x <- t(expr5)
y <- factor(label$Group)


set.seed(1012)

trainIndex <- createDataPartition(
  y,
  p = 0.7,
  list = FALSE
)

x_train <- x[trainIndex, ]
x_test  <- x[-trainIndex, ]

y_train <- y[trainIndex]
y_test  <- y[-trainIndex]


cvfit <- cv.glmnet(
  x_train,
  y_train,
  family = "binomial"
)


pdf("LASSO_CV_and_CoefficientPaths.pdf",
    width = 10, height = 5)

par(mfrow = c(1, 2))

plot(cvfit,
     xlab = "Log(lambda)",
     ylab = "Mean-Squared Error")

plot(cvfit$glmnet.fit,
     xvar = "lambda",
     label = TRUE,
     xlab = "Log(lambda)",
     ylab = "Coefficients")

dev.off()

cat("Save Lasso PDF：LASSO_CV_and_CoefficientPaths.pdf\n")


lasso_prob <- predict(
  cvfit,
  newx = x_test,
  s = "lambda.min",
  type = "response"
)



lasso_prob <- as.numeric(as.vector(lasso_prob)) 


print(paste("LASSO prob sample:", head(lasso_prob)))

roc_lasso <- roc(y_test, lasso_prob)
auc_lasso <- round(auc(roc_lasso), 3)



svm_model <- svm(
  x_train,
  y_train,
  kernel = "linear",
  probability = TRUE
)

svm_pred <- predict(
  svm_model,
  x_test,
  probability = TRUE
)


svm_prob <- attr(svm_pred, "probabilities")[, 2]
svm_prob <- as.numeric(as.vector(svm_prob))

print(paste("SVM prob sample:", head(svm_prob)))

roc_svm <- roc(y_test, svm_prob)
auc_svm <- round(auc(roc_svm), 3)


rf_model <- randomForest(
  x_train,
  y_train,
  ntree = 500,
  importance = TRUE
)


accuracy_results <- data.frame(
  Num_Features = 1:5,
  Accuracy = numeric(5)
)

for (i in 1:5) {
  rf_tmp <- randomForest(
    x_train[, 1:i, drop = FALSE],
    y_train,
    ntree = 500
  )
  accuracy_results$Accuracy[i] <- mean(rf_tmp$err.rate[, 1])
}

p_rf <- ggplot(accuracy_results, aes(x = Num_Features, y = Accuracy)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  theme_minimal() +
  labs(x = "Number of Features", y = "Accuracy") +
  theme(text = element_text(size = 14))

ggsave("RandomForest_Accuracy_Curve.pdf",
       plot = p_rf,
       width = 6, height = 5)

cat("Save RF：RandomForest_Accuracy_Curve.pdf\n")


rf_prob <- predict(
  rf_model,
  x_test,
  type = "prob"
)[, 2]


rf_prob <- as.numeric(as.vector(rf_prob))

print(paste("RF prob sample:", head(rf_prob)))

roc_rf <- roc(y_test, rf_prob)
auc_rf <- round(auc(roc_rf), 3)

auc_table <- data.frame(
  Model = c("LASSO", "SVM", "Random Forest"),
  AUC   = c(auc_lasso, auc_svm, auc_rf)
)

print(auc_table)


png(
  "ROC_5genes.png",
  width  = 3000,
  height = 2500,
  res     = 300
)

plot(
  roc_lasso,
  col = "#1B9E77",
  lwd = 3,
  legacy.axes = TRUE,
  main = "ROC Curves (5 Genes)"
)

lines(roc_svm, col = "#D95F02", lwd = 3)
lines(roc_rf,  col = "#7570B3", lwd = 3)

abline(a = 0, b = 1, lty = 2, col = "grey")

legend(
  "bottomright",
  legend = c(
    paste0("LASSO (AUC = ", auc_lasso, ")"),
    paste0("SVM (AUC = ", auc_svm, ")"),
    paste0("Random Forest (AUC = ", auc_rf, ")")
  ),
  col  = c("#1B9E77", "#D95F02", "#7570B3"),
  lwd  = 3,
  bty  = "n"
)

dev.off()
library(pROC)


ci_lasso <- ci.auc(roc_lasso, method = "delong")
ci_svm   <- ci.auc(roc_svm,   method = "delong")
ci_rf     <- ci.auc(roc_rf,     method = "delong")


print(ci_lasso)
print(ci_svm)
print(ci_rf)


png(
  "ROC_5genes.png",
  width  = 3000,
  height = 2500,
  res     = 300
)


par(family = "Times New Roman")

plot(
  roc_lasso,
  col = "#1B9E77",
  lwd = 3,
  legacy.axes = TRUE,
  main = "ROC Curves (5 Genes, GSE43292)",
  

  cex.main = 1.5,        
  font.main = 2,           
  

  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  cex.lab = 1.4,          
  font.lab = 2,           
  

  cex.axis = 1.2          
)

lines(roc_svm, col = "#D95F02", lwd = 3)
lines(roc_rf,  col = "#7570B3", lwd = 3)

abline(a = 0, b = 1, lty = 2, col = "grey", lwd = 2)

legend(
  "bottomright",
  legend = c(
    paste0("LASSO (AUC = ", auc_lasso, ", 95% CI ",
           sprintf("%.3f–%.3f", ci_lasso[1], ci_lasso[3]), ")"),
    paste0("SVM (AUC = ", auc_svm, ", 95% CI ",
           sprintf("%.3f–%.3f", ci_svm[1], ci_svm[3]), ")"),
    paste0("Random Forest (AUC = ", auc_rf, ", 95% CI ",
           sprintf("%.3f–%.3f", ci_rf[1], ci_rf[3]), ")")
  ),
  col  = c("#1B9E77", "#D95F02", "#7570B3"),
  lwd  = 3,
  bty  = "n",
  

  cex = 1.1,            
  text.font = 2            
)

dev.off()


coef_min <- coef(cvfit, s = "lambda.min")
coef_df <- as.data.frame(as.matrix(coef_min))


coef_df$Gene <- rownames(coef_df)


coef_df <- coef_df[coef_df$Gene != "(Intercept)", ]


coef_df <- coef_df[coef_df$lambda.min != 0, ]


colnames(coef_df)[1] <- "Coef"

print(coef_df)
write.csv(coef_df, "LASSO_selected_genes.csv", row.names = FALSE)

lasso_genes <- coef_df$Gene
w <- t(svm_model$coefs) %*% svm_model$SV

svm_imp <- data.frame(
  Gene = colnames(x_train),
  Weight = as.numeric(w)
)

svm_imp <- svm_imp[svm_imp$Gene %in% lasso_genes, ]
svm_imp <- svm_imp[order(-abs(svm_imp$Weight)), ]

print(svm_imp)

imp <- importance(rf_model)

imp_df <- data.frame(
  Gene = rownames(imp),
  Importance = imp[, 1]
)

imp_df <- imp_df[imp_df$Gene %in% lasso_genes, ]
imp_df <- imp_df[order(-imp_df$Importance), ]

print(imp_df)

library(VennDiagram)
library(readr)
library(tidyverse)


lasso_top3 <- coef_df %>%
  arrange(desc(abs(Coef))) %>%
  slice_head(n = 3) %>%
  pull(Gene)

print(lasso_top3)

write.table(
  lasso_top3,
  file = "LASSO_top3_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



svm_top3 <- svm_imp %>%
  arrange(desc(abs(Weight))) %>%
  slice_head(n = 3) %>%
  pull(Gene)

print(svm_top3)

write.table(
  svm_top3,
  file = "SVM_top3_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



rf_top3 <- imp_df %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 3) %>%
  pull(Gene)

print(rf_top3)

write.table(
  rf_top3,
  file = "RF_top3_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



library(VennDiagram)


ls_intersect <- intersect(lasso_top3, svm_top3)
print(ls_intersect)
write.table(
  ls_intersect,
  file = "Intersection_LASSO_SVM.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


venn.plot <- venn.diagram(
  x = list(
    LASSO = lasso_top3,
    SVM   = svm_top3
  ),
  filename = "Venn_LASSO_SVM_top3.tiff",
  
  fill = c("#E5C1B3", "#F1D27B"),
  alpha = 0.5,
  scaled = FALSE,
  euler.d = FALSE,
  
  lwd = c(0, 0),
  
  cex = 1.4,
  fontface = "bold",
  fontfamily = "Times New Roman",
  
  cat.cex = 1.2,
  cat.fontfamily = "Times New Roman",
  cat.fontface = "bold",
  cat.col = c("#C76A49", "#F1A700"),
  
  cat.pos = c(-6, 6),
  cat.dist = c(0.09, 0.09),
  cat.default.pos = "outer",
  
  margin = 0.20,
  
  insets = list(
    list(text = lr_intersect, cex = 0.85, fontface = "bold", col = "black")
  )
)

# LASSO ∩ Random Forest
lr_intersect <- intersect(lasso_top3, rf_top3)
print(lr_intersect)
write.table(
  lr_intersect,
  file = "Intersection_LASSO_RF.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


venn.diagram(
  x = list(LASSO = lasso_top3, RandomForest = rf_top3),
  filename = "Venn_LASSO_RF.tiff",
  fill = c("#E5C1B3", "#F3C6D3"),
  alpha = 0.5,
  scaled = FALSE,
  euler.d = FALSE,
  lwd = c(0, 0),
  cex = 1.4,
  fontface = "bold",
  fontfamily = "Times New Roman",
  cat.cex = 1.2,
  cat.fontfamily = "Times New Roman",
  cat.fontface = "bold",
  cat.col = c("#C76A49", "#C05A83"),
  cat.pos = c(-6, 6),
  cat.dist = c(0.09, 0.09),
  cat.default.pos = "outer",
  margin = 0.20
)


# SVM ∩ Random Forest
sr_intersect <- intersect(svm_top3, rf_top3)

write.table(
  sr_intersect,
  file = "Intersection_SVM_RF.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

print(sr_intersect)
venn.diagram(
  x = list(SVM = svm_top3, RandomForest = rf_top3),
  filename = "Venn_SVM_RF.tiff",
  fill = c("#F1D27B", "#F3C6D3"),
  alpha = 0.5,
  scaled = FALSE,
  euler.d = FALSE,
  lwd = c(0, 0),
  cex = 1.4,
  fontface = "bold",
  fontfamily = "Times New Roman",
  cat.cex = 1.2,
  cat.fontfamily = "Times New Roman",
  cat.fontface = "bold",
  cat.col = c("#F1A700", "#C05A83"),
  cat.pos = c(-6, 6),
  cat.dist = c(0.09, 0.09),
  cat.default.pos = "outer",
  margin = 0.20
)


