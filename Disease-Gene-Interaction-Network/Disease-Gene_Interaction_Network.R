library(tidyverse)
library(igraph)
library(ggraph)
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
    size = 3.5,
    box.padding = 0.6,   
    point.padding = 0.8, 
    force = 2,           
    max.overlaps = Inf
  ) +
  scale_size_manual(
    values = c(
      "Method" = 14,
      "Lasso+GBM" = 8,
      "MCODE" = 8
    )
  ) +
  scale_color_manual(
    values = c(
      "Method" = "#B8D8B8",
      "Lasso+GBM" = "#E5C1B3",
      "MCODE" = "#F1D27B"
    )
  ) +
  theme_void()


