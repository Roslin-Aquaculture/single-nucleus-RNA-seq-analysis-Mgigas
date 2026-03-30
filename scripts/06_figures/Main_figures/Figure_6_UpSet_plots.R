# Figure 6: differential expression outputs, unified heatmaps, and UpSet plots

# This script:
#   1. prepares the Seurat object for cluster- and stage-specific comparisons
#   2. runs differential expression analysis across infection stages
#   3. generates per-comparison plots and heatmaps
#   4. builds unified DE gene lists across stages
#   5. creates UpSet plots for genes upregulated in control or infection


############################### Step 1: DATA PREPARATION PRIOR TO DGE ##############################

############################
# Load required packages
############################

library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(scCustomize)
library(dittoSeq)
library(Seurat.utils)
library(ComplexUpset)
library(readr)

############################
# Load Seurat object and prepare metadata for group comparisons
############################

# Set working directory
setwd("/path/to/processed/seurat/object")

# Load processed Seurat object; this was created previously under scripts/05_seurat_processing/01_process_seurat_objects.R
seu_obj <- read_rds("seu_obj_filt_umap_18d_6r_3kRes_2026.rds") 

seu_obj

# Rename clusters using analysis-friendly labels with underscores
new.cluster.ids <- c(
  "Cluster_0", "Cluster_1",
  "Gill_ciliary_cells", "Hepatopancreas_cells",
  "Gill_neuroepithelial_cells", "Gill_cell_type_1",
  "Hyalinocytes", "Haemocyte_cell_type_1",
  "Mantle_cell_type_1", "Cluster_9",
  "Vesicular_haemocytes", "Immature_haemocytes",
  "Macrophage_like_cells", "Adductor_muscle_cells",
  "Mantle_cell_type_2", "Mantle_epithelial_cells",
  "Gill_cell_type_2", "Small_granule_cells"
)

names(new.cluster.ids) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new.cluster.ids)

############################
# Create simplified infection-stage groups and remove 24-hpiA sample
############################

# Add simplified infection-stage labels for downstream comparisons
seu_obj$condition_new <- case_when(
  seu_obj$sample %in% c("Homogenate", "Uninfected") ~ "control",
  seu_obj$sample %in% c("6-hpiA", "6-hpiD") ~ "6hpi",
  seu_obj$sample == "24-hpiA" ~ "mid?",
  seu_obj$sample == "24-hpiJ" ~ "24hpi",
  seu_obj$sample == "72-hpiJ" ~ "72hpi",
  seu_obj$sample == "96-hpiE" ~ "96hpi"
)

# Set factor order for downstream plotting and aggregation
seu_obj$condition_new <- factor(
  seu_obj$condition_new,
  levels = c("control", "6hpi", "mid?", "24hpi", "72hpi", "96hpi")
)

# Remove 24-hpiA ("mid?") from downstream analysis
seu_obj_clean <- subset(seu_obj, subset = condition_new != "mid?")
rm(seu_obj)

# Create a combined identity column containing cluster name and condition
# This is used for pairwise DE testing within each cluster across stages
seu_obj_clean$sample_subtypes <- paste(Idents(seu_obj_clean), seu_obj_clean$condition_new, sep = "_")
Idents(seu_obj_clean) <- "sample_subtypes"
unique(seu_obj_clean$sample_subtypes)

############################
# Add readable gene symbols to the Seurat object
############################

# Replace gene IDs with combined gene symbol / description labels
# to improve readability in downstream tables and plots

# Load annotation tables; these can be found in 06_figures/supporting_files
ORSON <- read_csv("01_ORSON_french_group_2024_biorxiv_suppl2.csv") %>%
  select(Gene.ID, Description = Sequence.Description)

cg_science <- read_tsv("02_Cg_gene_names.tsv", col_names = FALSE) %>%
  rename(Gene.ID = X1, Gene.symbol.science = X2)

# Merge annotations and build combined display labels
gene_map2 <- full_join(ORSON, cg_science, by = "Gene.ID") %>%
  mutate(
    gene_symbol = paste(
      coalesce(Gene.symbol.science, ""),
      coalesce(Description, ""),
      sep = " : "
    ) %>%
      trimws() %>%
      gsub("[^[:alnum:] :]", " ", x = .)
  ) %>%
  select(gene_id = Gene.ID, gene_symbol)

# Make duplicated gene labels unique
gene_map2 <- gene_map2 %>%
  group_by(gene_symbol) %>%
  mutate(
    gene_symbol = if (n() > 1) paste0(gene_symbol, LETTERS[seq_along(gene_symbol)]) else gene_symbol
  ) %>%
  ungroup()

# Prepend gene IDs to the display labels
gene_map2 <- gene_map2 %>%
  mutate(gene_symbol = paste(gene_id, gene_symbol, sep = " "))

# Keep only genes present in the Seurat object
gene_map_filtered <- gene_map2 %>%
  filter(gene_id %in% rownames(seu_obj_clean))

# Build rename vector for Seurat
gene_symbols <- gene_map_filtered$gene_symbol[match(rownames(seu_obj_clean), gene_map_filtered$gene_id)]
rename_vector <- setNames(gene_symbols, rownames(seu_obj_clean))
rename_vector[is.na(rename_vector)] <- names(rename_vector)[is.na(rename_vector)]
rename_vector <- make.unique(rename_vector)

# Rename genes in the Seurat object
seu_obj_clean <- RenameGenesSeurat(
  obj = seu_obj_clean,
  newnames = rename_vector
)

rownames(seu_obj_clean)

############################### DATA PREPARATION PRIOR TO DGE ENDS ##############################





###################################### STEP 2: DGE ANALYSIS #####################################

# Main output directory for all DE results and plots
parent_dir <- "de_plots_full_ann_30032026"
dir.create(parent_dir, showWarnings = FALSE)

# Set up log file
log_file <- file.path(parent_dir, "dge_analysis_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=====================================================================\n")
cat("Differential Expression Analysis Log\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("Parent directory:", parent_dir, "\n")
cat("=====================================================================\n\n")

# Define colour palette for heatmaps
pd_palette <- colorRampPalette(c(
  "#313695", "#4575B4", "#74ADD1",
  "#ABD9E9", "#E0F3F8", "#F7F7F7",
  "#FEE090", "#FDAE61", "#F46D43",
  "#D73027", "#A50026"
))

perform_comparison <- function(condition) {
  condition_dir <- file.path(parent_dir, condition)
  dir.create(condition_dir, showWarnings = FALSE)
  
  cat("\n=====================================================================\n")
  cat("Starting condition:", condition, "\n")
  cat("=====================================================================\n")
  
  # Identify control and target identities for the current condition
  control_clusters <- unique(seu_obj_clean$sample_subtypes[grepl("_control$", seu_obj_clean$sample_subtypes)])
  target_clusters  <- unique(seu_obj_clean$sample_subtypes[grepl(paste0("_", condition, "$"), seu_obj_clean$sample_subtypes)])
  
  control_tissues <- gsub("_control$", "", control_clusters)
  target_tissues  <- gsub(paste0("_", condition, "$"), "", target_clusters)
  
  for (i in seq_along(control_clusters)) {
    control_label <- control_clusters[i]
    tissue_name   <- control_tissues[i]
    target_label  <- paste0(tissue_name, "_", condition)
    
    cat("\n------------------------------\n")
    cat("Processing tissue:", tissue_name, "\n")
    cat("Control:", control_label, "\n")
    cat("Target:", target_label, "\n")
    cat("------------------------------\n")
    
    # Skip comparison if the matching target identity is absent
    if (!target_label %in% target_clusters) {
      cat("SKIP:", tissue_name, "- no matching target for", condition, "\n")
      next
    }
    
    comparison_dir <- file.path(condition_dir, paste0(tissue_name, "_control_vs_", condition))
    dir.create(comparison_dir, showWarnings = FALSE)
    
    # Run differential expression
    cli::cat_rule(paste("DGE:", control_label, "vs", target_label))
    
    cluster_control_target <- FindMarkers(
      seu_obj_clean,
      logfc.threshold = 1,
      min.pct = 0.25,
      ident.1 = control_label,
      ident.2 = target_label
    )
    
    # Handle comparisons with no DE output
    if (nrow(cluster_control_target) == 0) {
      cat("❌ No DE results for", tissue_name, "(", control_label, "vs", target_label, ")\n")
      next
    }
    
    cluster_control_target <- cluster_control_target %>%
      filter(p_val_adj < 0.01) %>%
      rownames_to_column(var = "gene")
    
    # Split DE genes into control-up and target-up sets
    up_in_control <- cluster_control_target %>% filter(avg_log2FC > 1) %>% pull(gene)
    up_in_target  <- cluster_control_target %>% filter(avg_log2FC < -1) %>% pull(gene)
    
    up_in_control <- up_in_control[!is.na(up_in_control) & up_in_control != ""]
    up_in_target  <- up_in_target[!is.na(up_in_target)  & up_in_target  != ""]
    
    cat(sprintf("✅ Found %d total DE genes (adj p < 0.01)\n", nrow(cluster_control_target)))
    cat(sprintf("   ↳ %d up in control | %d up in %s\n", length(up_in_control), length(up_in_target), condition))
    
    # Helper function to generate violin or dot plots for top DE genes
    plot_genes <- function(genes, plot_type, suffix, label_for_title) {
      if (length(genes) > 2) {
        genes_in_seurat <- genes[genes %in% rownames(seu_obj_clean)]
        if (length(genes_in_seurat) == 0) {
          cat("⚠️  No matching genes in Seurat object for", suffix, "\n")
          return(NULL)
        }
        
        cat("Plotting", plot_type, "for", tissue_name, "-", suffix, "(", length(genes_in_seurat), "genes)\n")
        fn <- file.path(comparison_dir, paste0(tissue_name, "_control_vs_", condition, "_", suffix, ".pdf"))
        pdf(fn, width = 12, height = 10)
        plot_title <- paste0(label_for_title, ": ", tissue_name, " (control vs ", condition, ")")
        if (plot_type == "VlnPlot") {
          print(
            VlnPlot(
              seu_obj_clean,
              features = genes_in_seurat[1:min(10, length(genes_in_seurat))],
              stack = TRUE,
              flip = TRUE,
              group.by = "condition_new"
            ) +
              ggtitle(plot_title) +
              NoLegend()
          )
        } else {
          print(
            DotPlot(
              seu_obj_clean,
              features = rev(genes_in_seurat[1:min(10, length(genes_in_seurat))]),
              group.by = "condition_new"
            ) +
              coord_flip() +
              scale_x_discrete(position = "top") +
              ggtitle(plot_title)
          )
        }
        dev.off()
        cat("   → Saved plot:", fn, "\n")
      } else {
        cat("⚠️  Not enough genes for", suffix, "(need >2)\n")
      }
    }
    
    # Generate summary plots for DE genes
    plot_genes(up_in_control, "VlnPlot", "vln_up_in_control", "Top markers up in control")
    plot_genes(up_in_target,  "VlnPlot", "vln_up_in_target",  paste0("Top markers up in ", condition))
    plot_genes(up_in_control, "DotPlot", "dot_up_in_control", "Top markers up in control")
    plot_genes(up_in_target,  "DotPlot", "dot_up_in_target",  paste0("Top markers up in ", condition))
    
    # Aggregate expression across infection stages for heatmap plotting
    cluster_aggregate_exp <- AggregateExpression(seu_obj_clean, return.seurat = TRUE, group.by = "condition_new")
    cluster_aggregate_exp$condition_new <- factor(
      cluster_aggregate_exp$condition_new,
      levels = c("control", "6hpi", "24hpi", "72hpi", "96hpi")
    )
    
    expr_mat <- GetAssayData(cluster_aggregate_exp, slot = "data")
    expr_mat_subset <- expr_mat[cluster_control_target$gene, , drop = FALSE]
    expr_df <- as.data.frame(as.matrix(expr_mat_subset))
    expr_df <- tibble::rownames_to_column(expr_df, var = "gene")
    comparison_label <- paste0(tissue_name, "_vs_", condition)
    colnames(expr_df)[-1] <- paste0(comparison_label, ".", colnames(expr_df)[-1])
    expr_path <- file.path(comparison_dir, paste0(comparison_label, "_heatmap_expr_table.tsv"))
    write_tsv(expr_df, expr_path)
    cat("Saved expression table:", expr_path, "\n")
    
    # Generate DE heatmap
    heatmap_file <- file.path(comparison_dir, paste0(tissue_name, "_control_vs_", condition, "_de_genes_heatmap.pdf"))
    cat("Generating heatmap for", tissue_name, "(", condition, ")\n")
    heatmap_plot <- dittoHeatmap(
      cluster_aggregate_exp,
      genes = cluster_control_target$gene,
      annot.by = "condition_new",
      heatmap.colors = pd_palette(256),
      main = paste0(
        "Total ", length(cluster_control_target$gene), " DE genes: ",
        tissue_name, " (control vs ", condition, ")"
      ),
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      scale = "row",
      show_colnames = TRUE
    )
    pdf(heatmap_file, width = 12, height = 10)
    print(heatmap_plot)
    dev.off()
    cat("   → Saved heatmap:", heatmap_file, "\n")
    
    # Save DE gene lists
    base_filename <- paste0(tissue_name, "_control_vs_", condition, "_de_genes_")
    write.table(
      cluster_control_target$gene,
      file.path(comparison_dir, paste0(base_filename, "all.txt")),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
    write.table(
      up_in_control,
      file.path(comparison_dir, paste0(base_filename, "up_in_control.txt")),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
    write.table(
      up_in_target,
      file.path(comparison_dir, paste0(base_filename, "up_in_target.txt")),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
    
    cat("Saved DE gene lists to:", comparison_dir, "\n")
  }
}

# Run pairwise comparisons for all infection stages
comparison_conditions <- c("6hpi", "24hpi", "72hpi", "96hpi")
lapply(comparison_conditions, perform_comparison)

cat("End time:", as.character(Sys.time()), "\n")
cat("Log saved to:", log_file, "\n")
cat("=====================================================================\n")

sink(type = "message")
sink(type = "output")
close(log_con)

###################################### DGE ANALYSIS ENDS #################################



######################## STEP 3: UNIFIED GENE LISTS AND HEATMAPS ########################

# Create output directories for unified gene lists and heatmaps
unified_dir <- file.path(parent_dir, "unified_gene_lists")
dir.create(unified_dir, showWarnings = TRUE)

heatmap_dir <- file.path(unified_dir, "heatmap_plots")
dir.create(heatmap_dir, showWarnings = FALSE)

# Set up log file for unified heatmap generation
log_file <- file.path(unified_dir, "unified_heatmap_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("=====================================================================\n")
cat("Unified Gene List & Heatmap Generation Log\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("=====================================================================\n\n")

# Aggregate expression once for all unified heatmap plots
cluster_aggregate_exp <- AggregateExpression(seu_obj_clean, return.seurat = TRUE, group.by = "condition_new")

# Set infection-stage order
cluster_aggregate_exp$condition_new <- factor(
  cluster_aggregate_exp$condition_new,
  levels = c("control", "6hpi", "24hpi", "72hpi", "96hpi")
)

expr_mat <- GetAssayData(cluster_aggregate_exp, slot = "data")

############################################################
# Locate all DE gene list files
de_files <- list.files(parent_dir, pattern = "_de_genes_all.txt$", recursive = TRUE, full.names = TRUE)

# Initialize list to store DE genes by cluster comparison
cluster_gene_list <- list()

# Populate the list of DE genes for each cluster comparison
for (file in de_files) {
  cluster_name <- basename(dirname(file))
  genes <- readLines(file)
  
  if (cluster_name %in% names(cluster_gene_list)) {
    cluster_gene_list[[cluster_name]] <- unique(c(cluster_gene_list[[cluster_name]], genes))
  } else {
    cluster_gene_list[[cluster_name]] <- unique(genes)
  }
}

################################################################
# Derive base cluster names prior to the "_control_vs_*" suffix
processed_clusters <- unique(sub("_control_vs_.*", "", names(cluster_gene_list)))

for (base_name in processed_clusters) {
  base_name_clean <- base_name
  combined_genes <- character()
  
  cat("\n==============================\n")
  cat("Processing cluster:", base_name_clean, "\n")
  cat("==============================\n")
  
  for (condition in c("6hpi", "24hpi", "72hpi", "96hpi")) {
    full_cluster_name <- paste0(base_name, "_control_vs_", condition)
    
    if (full_cluster_name %in% names(cluster_gene_list)) {
      n_genes <- length(cluster_gene_list[[full_cluster_name]])
      cat("  ✅ Found:", full_cluster_name, " | ", n_genes, " genes\n")
      combined_genes <- unique(c(combined_genes, cluster_gene_list[[full_cluster_name]]))
    } else {
      cat("  ⚠️ Missing:", full_cluster_name, "\n")
    }
  }
  
  cat("  → Total combined genes for", base_name_clean, ":", length(combined_genes), "\n")
  
  # Skip if no DE gene files were found for the cluster
  if (length(combined_genes) == 0) {
    warning(paste("No DE genes found for cluster:", base_name_clean))
    next
  }
  
  # Save unified gene list
  write.table(
    combined_genes,
    file = file.path(unified_dir, paste0(base_name_clean, "_unified_genes.txt")),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  # Keep only genes present in the aggregated expression matrix
  valid_genes <- combined_genes[combined_genes %in% rownames(expr_mat)]
  
  cat("  → Valid genes present in expression matrix:", length(valid_genes), "\n")
  
  if (length(valid_genes) == 0) {
    warning(paste("No valid genes found in expression matrix for:", base_name_clean))
    next
  }
  
  cat("  🔥 Generating heatmap for", base_name_clean, "...\n")
  
  # Generate unified heatmap
  heatmap_plot <- dittoHeatmap(
    object = cluster_aggregate_exp,
    genes = valid_genes,
    annot.by = "condition_new",
    order.by = "condition_new",
    heatmap.colors = pd_palette(256),
    main = paste0("Unified DE genes heatmap: ", base_name_clean, " (", length(valid_genes), " genes)"),
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    scale = "row",
    show_colnames = TRUE
  )
  
  # Save unified heatmap as PDF
  heatmap_file <- paste0(base_name_clean, "_unified_genes_heatmap.pdf")
  pdf(file.path(heatmap_dir, heatmap_file), width = 12, height = 10)
  print(heatmap_plot)
  dev.off()
  
  cat("  ✅ Saved heatmap to:", file.path(heatmap_dir, heatmap_file), "\n")
}

cat("\n=====================================================================\n")
cat("✅ All unified gene lists and heatmaps generated successfully.\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("Log saved to:", log_file, "\n")
cat("=====================================================================\n")

sink(type = "message")
sink(type = "output")
close(log_con)

######################## UNIFIED GENE LISTS AND HEATMAPS END ################################




######################## STEP 4: CREATE UPSET PLOTS #######################################

#............. Genes upregulated in control: STARTS ........................................

stages <- c("6hpi", "24hpi", "72hpi", "96hpi")

# Load all control-upregulated gene lists and keep stage as metadata
combined_up_all <- map_dfr(stages, function(stage) {
  files <- list.files(
    path = file.path(parent_dir, stage),
    pattern = "_de_genes_up_in_control.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(files) == 0) return(NULL)
  
  map_dfr(files, function(f) {
    cluster <- basename(dirname(f))
    read_tsv(f, col_names = "gene", show_col_types = FALSE) %>%
      mutate(cluster = cluster, stage = stage)
  })
})

# Remove comparison-specific suffixes from cluster labels
combined_up_all_clean <- combined_up_all %>%
  mutate(
    cluster = str_remove(cluster, "_control_vs_.*$")
  )

# Create logical matrix for ComplexUpset
gene_matrix_up <- combined_up_all_clean %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    names_from = cluster,
    values_from = present,
    values_fill = FALSE
  ) %>%
  select(gene, stage, everything())

# Set infection-stage order
gene_matrix_up <- gene_matrix_up %>%
  mutate(stage = factor(stage, levels = c("6hpi", "24hpi", "72hpi", "96hpi")))

clusters <- setdiff(colnames(gene_matrix_up), c("gene", "stage"))

# Save UpSet plot for genes upregulated in control
pdf(
  file = file.path(parent_dir, "Upregulated_in_Control_UpSetPlot.pdf"),
  width = 20,
  height = 14
)

stage_colors <- c(
  `6hpi`   = "#1b9e77",
  `24hpi`  = "#d95f02",
  `72hpi`  = "#7570b3",
  `96hpi`  = "#386cb0"
)

upset(
  gene_matrix_up,
  intersect = clusters,
  name = "Upregulated in control",
  min_size = 10,
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = stage)
    ) +
      scale_fill_manual(values = stage_colors) +
      labs(fill = "Stage") +
      theme_minimal(base_size = 22) +
      theme(
        legend.position = "top",
        legend.title = element_text(face = "bold")
      )
  )
)

dev.off()

#............... Genes upregulated in control: ENDS .............................




#............. Genes upregulated in infected samples: STARTS ....................

# Load all infection-upregulated gene lists and keep stage as metadata
combined_down_all <- map_dfr(stages, function(stage) {
  files <- list.files(
    path = file.path(parent_dir, stage),
    pattern = "_de_genes_up_in_target.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(files) == 0) return(NULL)
  
  map_dfr(files, function(f) {
    cluster <- basename(dirname(f))
    read_tsv(f, col_names = "gene", show_col_types = FALSE) %>%
      mutate(cluster = cluster, stage = stage)
  })
})

# Remove comparison-specific suffixes from cluster labels
combined_down_all_clean <- combined_down_all %>%
  mutate(
    cluster = str_remove(cluster, "_control_vs_.*$")
  )

# Create logical matrix for ComplexUpset
gene_matrix_down <- combined_down_all_clean %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    names_from = cluster,
    values_from = present,
    values_fill = FALSE
  ) %>%
  select(gene, stage, everything())

# Set infection-stage order
gene_matrix_down <- gene_matrix_down %>%
  mutate(stage = factor(stage, levels = c("6hpi", "24hpi", "72hpi", "96hpi")))

clusters <- setdiff(colnames(gene_matrix_down), c("gene", "stage"))

# Save UpSet plot for genes upregulated in infection
pdf(
  file = file.path(parent_dir, "Upregulated_in_Infection_UpSetPlot.pdf"),
  width = 20,
  height = 14
)

upset(
  gene_matrix_down,
  intersect = clusters,
  name = "Upregulated in infection",
  min_size = 10,
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = stage)
    ) +
      scale_fill_manual(values = stage_colors) +
      labs(fill = "Stage") +
      theme_minimal(base_size = 22) +
      theme(
        legend.position = "top",
        legend.title = element_text(face = "bold")
      )
  )
)

dev.off()

#............. Genes upregulated in infected samples: ENDS .....................

# ------------------------------------------------------------------
# End of Figure 6 script
# ------------------------------------------------------------------
