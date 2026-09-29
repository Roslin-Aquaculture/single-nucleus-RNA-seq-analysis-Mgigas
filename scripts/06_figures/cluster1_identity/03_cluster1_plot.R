# =========================================================
# 07d_cluster1_identity_figure.R
# Figure: terms enriched in Cluster 1 and in no other cluster
# Environment: the one with ggplot2 / svglite. No clusterProfiler needed.
# =========================================================
# ASCII only throughout.
#
# Recomputes nothing. Reads cluster1_discriminating_terms.tsv, written by
# 07b_cluster1_identity_enrichment.R, so the two scripts can run in
# different environments.
#
# The figure is the visual form of the claim in the identity paragraph:
# terms enriched in Cluster 1's control-only markers that are enriched in
# no other cluster's.
# =========================================================

library(readr); library(dplyr)
library(ggplot2); library(svglite); library(systemfonts); library(xml2)

# =========================================================
# CONFIG
# =========================================================
IN_DIR <- "/home/pdewari/eggnog/results/plots/cluster1_identity_280926"
IN_TSV <- file.path(IN_DIR, "cluster1_discriminating_terms.tsv")
STEM   <- file.path(IN_DIR, "FigS_cluster1_discriminating_terms")

N_SHOW      <- 25       # terms drawn, most significant first
EXCLUSIVE   <- TRUE     # TRUE: only terms shared with no other cluster

FONT       <- "Arial"   # Liberation Sans substitutes on Linux
BASE_SIZE  <- 8
PANEL_W_MM <- 150
PANEL_H_MM <- 180
WRAP_AT    <- 40        # characters before a term label wraps

fm <- systemfonts::match_fonts(FONT)
cat("font requested:", FONT, "-> resolved:", basename(fm$path), "\n")

# =========================================================
# EXPORT HELPERS (same standard as 07 / 08 / the IAP figure)
# =========================================================
# svglite puts shared stroke and fill declarations in a CSS block inside
# <defs>; Inkscape drops that stylesheet on ungroup and axis lines, ticks
# and panel borders disappear. It also writes textLength on every <text>,
# so deleting a character in Inkscape stretches the rest to fill the
# original width. Both are undone here.
fix_svg_for_inkscape <- function(file) {
  x  <- xml2::read_xml(file); ns <- xml2::xml_ns(x)
  els <- xml2::xml_find_all(
    x, "//d1:line|//d1:polyline|//d1:polygon|//d1:path|//d1:rect|//d1:circle", ns)
  n <- 0
  for (e in els) {
    s <- xml2::xml_attr(e, "style"); if (is.na(s)) s <- ""
    changed <- FALSE
    if (!grepl("stroke\\s*:", s)) {
      s <- paste0(trimws(s), if (nzchar(trimws(s))) " " else "", "stroke: #000000;")
      changed <- TRUE
    }
    if (!grepl("fill\\s*:", s)) { s <- paste0(trimws(s), " fill: none;"); changed <- TRUE }
    if (changed) { xml2::xml_set_attr(e, "style", s); n <- n + 1 }
  }
  tx <- xml2::xml_find_all(x, "//d1:text|//d1:tspan", ns)
  for (e in tx) {
    xml2::xml_set_attr(e, "textLength",   NULL)
    xml2::xml_set_attr(e, "lengthAdjust", NULL)
  }
  xml2::write_xml(x, file)
  cat("     inlined styles on", n, "of", length(els), "elements;",
      "released", length(tx), "text widths\n")
}

# cairo_pdf keeps text as text and carries physical units unambiguously.
# The default pdf() device writes per-glyph positions and splits letters
# apart in Inkscape.
save_both <- function(plot, stem, width_mm, height_mm) {
  ggsave(paste0(stem, ".svg"), plot, width = width_mm, height = height_mm,
         units = "mm", device = svglite::svglite, limitsize = FALSE)
  fix_svg_for_inkscape(paste0(stem, ".svg"))
  ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, limitsize = FALSE)
}

wrap_lab <- function(x, n = WRAP_AT) {
  vapply(x, function(s) paste(strwrap(s, n), collapse = "\n"), character(1))
}

# =========================================================
# LOAD
# =========================================================
stopifnot(file.exists(IN_TSV))
disc <- read_tsv(IN_TSV, show_col_types = FALSE)

cat("terms read:", nrow(disc),
    "| shared with no other cluster:", sum(disc$n_other_clusters_sharing == 0), "\n")

top_disc <- disc %>%
  { if (EXCLUSIVE) filter(., n_other_clusters_sharing == 0) else . } %>%
  arrange(p.adjust) %>%
  slice_head(n = N_SHOW) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

if (!nrow(top_disc)) stop("nothing to plot at these settings")

cat(sprintf("plotting %d terms (%g x %g mm; about %.0f mm high would fit them)\n",
            nrow(top_disc), PANEL_W_MM, PANEL_H_MM, 40 + 6 * nrow(top_disc)))

# =========================================================
# PLOT
# =========================================================
p_disc <- ggplot(top_disc,
                 aes(fold_enrichment, Description,
                     colour = -log10(p.adjust), size = Count)) +
  geom_point() +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = wrap_lab) +
  scale_colour_gradient(low = "#FCBBA1", high = "#A50F15",
                        limits = c(0, NA),
                        name   = "-log10 FDR") +
  scale_size_continuous(range = c(1.5, 5), name = "Genes") +
  guides(colour = guide_colourbar(frame.colour = NA, ticks = FALSE)) +
  labs(title = "Terms enriched in Cluster 1 and in no other cluster",
       x = "Fold enrichment", y = NULL) +
  theme_bw(base_size = BASE_SIZE, base_family = FONT) +
  theme(
    # every text element set explicitly; ggplot defaults several to grey30
    text              = element_text(size = BASE_SIZE, colour = "black", family = FONT),
    plot.title        = element_text(size = BASE_SIZE, colour = "black",
                                     family = FONT, face = "plain"),
    axis.text         = element_text(size = BASE_SIZE, colour = "black"),
    axis.title        = element_text(size = BASE_SIZE, colour = "black"),
    strip.text        = element_text(size = BASE_SIZE, colour = "black"),
    legend.text       = element_text(size = BASE_SIZE, colour = "black"),
    legend.title      = element_text(size = BASE_SIZE, colour = "black"),
    legend.key.height = unit(12, "pt"),
    legend.key.width  = unit(7,  "pt"),
    strip.background  = element_rect(fill = "grey95", colour = NA),
    panel.border      = element_rect(colour = "black", linewidth = 0.3, fill = NA),
    panel.grid.major  = element_line(colour = "grey94", linewidth = 0.25),
    panel.grid.minor  = element_blank(),
    axis.ticks        = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(1.5, "pt"),
    legend.position   = "right")

save_both(p_disc, STEM, PANEL_W_MM, PANEL_H_MM)

write_tsv(top_disc, paste0(STEM, "_plotted_terms.tsv"), na = "")
writeLines(capture.output(sessionInfo()), file.path(IN_DIR, "sessionInfo_07d.txt"))

cat("\nOutputs:", IN_DIR, "\n")
