install.packages("UpSetR")
library(UpSetR)
library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("D:\\R") 


gene_data <- read_excel("input.xlsx", sheet = "Sheet1")

# Create a gene list
gene_lists <- list(
  PFOS = gene_data$PFOS %>% na.omit() %>% as.character() %>% unique(),
  Stroke = gene_data$Stroke %>% na.omit() %>% as.character() %>% unique(),
  `Chronic kidney Disease` = gene_data$`Chronic kidney Disease` %>% na.omit() %>% as.character() %>% unique(),
  `Kidney Failure` = gene_data$`Kidney Failure` %>% na.omit() %>% as.character() %>% unique(),
  Hypertriglyceridemia = gene_data$Hypertriglyceridemia %>% na.omit() %>% as.character() %>% unique(),
  Hypertension = gene_data$Hypertension %>% na.omit() %>% as.character() %>% unique(),
  `Type Ⅱ Diabetes` = gene_data$`Type Ⅱ Diabetes` %>% na.omit() %>% as.character() %>% unique(),
  Obesity = gene_data$Obesity %>% na.omit() %>% as.character() %>% unique(),
  `Coronary Heart Disease` = gene_data$`Coronary Heart Disease` %>% na.omit() %>% as.character() %>% unique(),
  `Heart Failure` = gene_data$`Heart Failure` %>% na.omit() %>% as.character() %>% unique(),
  `Atherosclerotic Cardiovascular Disease` = gene_data$`Atherosclerotic Cardiovascular Disease` %>% na.omit() %>% as.character() %>% unique(),
  `Atrial Fibrillation` = gene_data$`Atrial Fibrillation` %>% na.omit() %>% as.character() %>% unique(),
  `Peripheral Artery Disease` = gene_data$`Peripheral Artery Disease` %>% na.omit() %>% as.character() %>% unique()
)

# Keep only the PFOS-intersecting genes
pfos_genes <- gene_lists$PFOS

# Filter out the other sets that share common genes with PFOS.
intersect_sets <- sapply(gene_lists, function(x) any(x %in% pfos_genes))
gene_lists_filtered <- gene_lists[intersect_sets]

# Filter the data to retain only the genes that intersect with PFOS.
filtered_genes <- unique(unlist(gene_lists))
filtered_genes <- filtered_genes[filtered_genes %in% pfos_genes]

# Create a new list for each set, containing only the genes that intersect with PFOS.
gene_lists_pfos <- lapply(gene_lists, function(gene_set) {
  intersect(gene_set, filtered_genes)
})


set_colors <- colorRampPalette(brewer.pal(12, "Set3"))(13)

# Plot an UpSet diagram
mat <- fromList(gene_lists_pfos)
library(dplyr)

# Combine the 0/1 values in each row into a key string.
comb_df <- mat %>%
  mutate(comb = apply(., 1, paste, collapse = "")) %>%
  group_by(comb) %>%
  summarise(freq = n()) %>%
  arrange(desc(freq))

print(comb_df)

target_comb <- comb_df$comb[comb_df$freq == 482]
print(target_comb)

# Extract the set names.
set_names <- colnames(mat)

# Convert the string to 0/1
bits <- strsplit(target_comb, "")[[1]] %>% as.numeric()

highlight_sets <- set_names[bits == 1]
print(highlight_sets)


#Generate the plot.

upset_plot <- upset(
  fromList(gene_lists_pfos),  # Use the PFOS-overlapping gene list
  nsets = 13,  
  nintersects = 30,  # Show top 30 intersections only
  mb.ratio = c(0.55, 0.45),
  order.by = "freq",
  decreasing = TRUE,
  text.scale = c(1.95, 1.95, 1.5, 1.5, 2.25, 1.8),  
  point.size = 4,
  line.size = 1,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Number of Genes per Set",
  set_size.show = TRUE,
  set_size.scale_max = max(sapply(gene_lists_pfos, length)) * 1.1,
  sets.bar.color = set_colors, 
  main.bar.color = "#87CEEB",  
  
  # Highlight intersections with PFOS
  queries = list(
    list(
      query = intersects,
      params = list(highlight_sets),
      color = "#E74C3C",
      active = TRUE
    )
  )
) 


print(upset_plot)
