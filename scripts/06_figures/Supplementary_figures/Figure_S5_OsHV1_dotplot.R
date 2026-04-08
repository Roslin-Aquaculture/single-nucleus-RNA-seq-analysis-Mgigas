# Figure S5: dotplot of viral transcript expression across infection stages and clusters

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(plotly)
library(shiny)


############################
# Load processed Seurat object
############################

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object generated previously in scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds")

############################
# Rename clusters with final cell-type annotations
############################

new.cluster.ids <- c(
  "Cluster 0", "Cluster 1",
  "Gill ciliary cells", "Hepatopancreas cells",
  "Gill neuroepithelial cells", "Gill cell type 1",
  "Hyalinocytes", "Haemocyte cell type 1",
  "Mantle cell type 1", "Cluster 9",
  "Vesicular haemocytes", "Immature haemocytes",
  "Macrophage-like cells", "Adductor muscle cells",
  "Mantle cell type 2", "Mantle epithelial cells",
  "Gill cell type 2", "Small granule cells"
)

# Match new labels to the current cluster levels and rename identities
names(new.cluster.ids) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new.cluster.ids)

############################
# Create simplified infection-stage groups
############################

# Collapse individual samples into broader infection-stage categories
seu_obj$condition_new <- dplyr::case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample == "24-hpiA" ~ "mid?",
  seu_obj$sample == "24-hpiJ" ~ "24hpi",
  seu_obj$sample == "72-hpiJ" ~ "72hpi",
  seu_obj$sample == "96-hpiE" ~ "96hpi"
)

# Set stage order before filtering
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove the 24-hpiA sample labelled as "mid?" from downstream plotting
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")
rm(seu_obj)

# Drop unused factor levels and keep the final plotting order
seu_obj_clean$condition_new <- factor(
  seu_obj_clean$condition_new,
  levels = c("control", "6hpi", "24hpi", "72hpi", "96hpi")
)

############################
# Define viral genes for plotting
############################

# Select viral ORFs from the feature names and exclude entries containing dots
gene_names <- rownames(seu_obj_clean)[
  grepl("^ORF", rownames(seu_obj_clean)) &
    !grepl("\\.", rownames(seu_obj_clean))
]

# Use ORF names directly as x-axis labels
alt_labels <- gene_names

############################
# Build combined condition-cluster grouping for dotplot display
############################

# Store cluster identities in metadata for downstream grouping
seu_obj_clean$cluster_name <- as.character(Idents(seu_obj_clean))

# Preserve the cluster order currently stored in the Seurat object
cluster_order <- levels(Idents(seu_obj_clean))

# Create a combined label so rows are arranged by condition, then cluster
seu_obj_clean$condition_cluster <- paste(
  seu_obj_clean$condition_new,
  seu_obj_clean$cluster_name,
  sep = "_"
)

# Build the desired display order:
# all clusters in control, then all clusters in 6hpi, then 24hpi, 72hpi, and 96hpi
combined_order <- unlist(lapply(levels(seu_obj_clean$condition_new), function(cond) {
  paste(cond, cluster_order, sep = "_")
}))

# Keep only combinations that are present in the object
combined_order <- combined_order[
  combined_order %in% unique(seu_obj_clean$condition_cluster)
]

# Apply the ordered factor levels
seu_obj_clean$condition_cluster <- factor(
  seu_obj_clean$condition_cluster,
  levels = combined_order
)

# Format row labels for plotting
clean_condition_cluster_labels <- sub(
  "^([^_]+)_(.*)$",
  "\\1 | \\2",
  levels(seu_obj_clean$condition_cluster)
)

############################
# Generate viral ORF dotplot
############################

p <- DotPlot(
  seu_obj_clean,
  features = gene_names,
  group.by = "condition_cluster",
  dot.scale = 8
) +
  scale_x_discrete(labels = alt_labels) +
  scale_y_discrete(labels = clean_condition_cluster_labels) +
  scale_color_gradientn(
    colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_blank(),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
  ) +
  coord_flip(clip = "off")

p
############################
# Save dotplot as PDF
############################

ggsave(
  filename = "Figure5_viral_ORF_dotplot.pdf",
  plot = p,
  width = 20,
  height = 14,
  device = cairo_pdf
)

############################
# Generate interactive plotly version
############################

ggplotly(p)


############################
# Generate interactive shiny version and search for an ORF
############################
ui <- fluidPage(
  titlePanel("Search viral ORF expression"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("orf", "Enter ORF name", value = "ORF100"),
      helpText("Example: ORF100, ORF112, ORF113"),
      width = 3
    ),
    
    mainPanel(
      plotlyOutput("orf_plot", height = "900px"),
      br(),
      verbatimTextOutput("msg")
    )
  )
)

server <- function(input, output, session) {
  
  selected_gene <- reactive({
    toupper(trimws(input$orf))
  })
  
  output$msg <- renderText({
    g <- selected_gene()
    
    if (!(g %in% gene_names)) {
      paste0(
        g, " not found. Check the ORF name in gene_names."
      )
    } else {
      paste0("Showing expression for ", g)
    }
  })
  
  output$orf_plot <- renderPlotly({
    g <- selected_gene()
    req(g %in% gene_names)
    
    p_gene <- DotPlot(
      seu_obj_clean,
      features = g,
      group.by = "condition_cluster",
      dot.scale = 8
    ) +
      scale_x_discrete(labels = g) +
      scale_y_discrete(labels = clean_condition_cluster_labels) +
      scale_color_gradientn(
        colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        text = element_text(family = "Arial"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title = element_blank(),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
      ) +
      coord_flip(clip = "off")
    
    ggplotly(p_gene)
  })
}

shinyApp(ui = ui, server = server)

############################
# Figure S5 ENDS 
############################
