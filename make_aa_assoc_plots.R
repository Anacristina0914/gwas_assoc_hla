library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

# Load input/output from snakemake
input_files = snakemake@input[["aa_labs"]]
output_summary = snakemake@output[["summary_table"]]
out_type1 = snakemake@output[["type1_aa_plot"]]
out_type2 = snakemake@output[["type2_aa_plot"]]
out_top_hm_t1 = snakemake@output[["hm_type1_plot"]]
out_top_hm_t2 = snakemake@output[["hm_type2_plot"]]

# Load groups colors
group_colors <- unlist(snakemake@params[["group_cols"]])

# Load aa dataframes
aa_all <- setNames(
  lapply(input_files, read_delim, delim = "\t", col_types = cols()),
  sub("\\.aa_labs\\.assoc.logistic$", "", basename(input_files))
)

sig_colors <- sort(group_colors[names(group_colors) != "nonsig"])

# Correct — sorts by the key names e.g. cluster1, cluster2, cluster3...
sig_color_names <- sort(names(group_colors[names(group_colors) != "nonsig"]))
cluster_names   <- sort(names(aa_all))

if (length(sig_color_names) != length(cluster_names)) {
  stop(sprintf(
    "Number of sig_* colors provided (%d) doesn't match number of cluster files (%d).\n  Expected: %s\n  Provided: %s",
    length(cluster_names),
    length(sig_color_names),
    paste(cluster_names, collapse = ", "),
    paste(sig_color_names, collapse = ", ")
  ))
}

if (!identical(sig_color_names, cluster_names)) {
  stop(sprintf(
    "Color name mismatch between config and cluster files.\n  Expected: %s\n  Provided: %s",
    paste(cluster_names, collapse = ", "),
    paste(sig_color_names, collapse = ", ")
  ))
}

# Add cluster name to each dataset, adjust pvalue,and add labels for class I and class II
aa_all <- lapply(names(aa_all), function(cl) {
  aa_all[[cl]] %>%
    mutate(
      cluster    = cl,
      adj_p      = p.adjust(P, method = "fdr"),
      log10Pval  = -log10(adj_p),
      significant = adj_p < 0.05,
      class      = ifelse(grepl("D[PQR][AB]1", locus), "ClassII", "ClassI")
    )
})

# Bind all datasets together
aa_all <- bind_rows(aa_all)

# Add cluster significant labels
aa_all <- aa_all %>%
  mutate(
    significant_cluster = ifelse(
      significant,
      cluster,
      "nonsig"
    ),
    variant = paste0(locus, ":", position, "_", aa)
  )

# Number of points to annotate
top_n <- snakemake@params[["top_n"]]

aa_all_top <- aa_all %>%
  group_by(cluster, locus) %>%
  mutate(
    rank_pval = rank(-log10Pval, ties.method = "min"),
    # Only label one point per position+pval combination to avoid stacked label duplicates
    is_top = significant & rank_pval <= top_n,
    dup_label = duplicated(paste(position, log10Pval)) & is_top,
    label = ifelse(is_top & !dup_label, residue, NA)
  ) %>%
  select(-rank_pval, -dup_label) %>%
  ungroup()


# ── Shared plot theme ─────────────────────────────────────────────────────────
xtop_size <- 16-(max(unlist(lapply(names(aa_all), nchar)))/1.7) # Set x-axis text size based on the longest var name

base_theme <- list(
  geom_point(alpha = 0.7, size = 1.5),
  geom_text_repel(
  aes(label = label),
  size               = 2.5,
  max.overlaps       = 20,
  na.rm              = TRUE,
  box.padding        = 0.4,
  point.padding      = 0.3,
  min.segment.length = 0,       # always draw the segment line
  segment.color      = "grey40",
  segment.size       = 0.3,
  show.legend        = FALSE
  ),
  scale_color_manual(values = group_colors),
  facet_grid(locus ~ cluster, scales = "free_x"),
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2),
  theme_bw(),
  theme(
    strip.background  = element_rect(fill = "white"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    strip.text.x.top  = element_text(face = "bold", size = xtop_size),
    #axis.title.y      = element_blank(),
    axis.text.y       = element_text(size = 8),
    axis.text.x       = element_text(size = 6),
    axis.title.x      = element_blank(),
    panel.spacing     = unit(0.8, "lines")
  )
)

# ── Class I plot ──────────────────────────────────────────────────────────────

p_type1 <- aa_all_top  %>%
  filter(class == "ClassI") %>%
  ggplot(aes(x = position, y = log10Pval, color = significant_cluster)) +
  base_theme +
  labs(
    title = "HLA Class I Associations",
    y     = "-log10(FDR corrected p-value)"
  ) + theme(legend.position="none")

ggsave(out_type1, plot = p_type1, device = "jpeg", width = 12, height = 6, dpi = 300)

# ── Class II plot ─────────────────────────────────────────────────────────────

p_type2 <- aa_all_top %>%
  filter(class == "ClassII") %>%
  ggplot(aes(x = position, y = log10Pval, color = significant_cluster)) +
  base_theme +
  labs(
    title = "HLA Class II Associations",
    y     = "-log10(FDR corrected p-value)"
  ) + theme(legend.position="none")

ggsave(out_type2, plot = p_type2, device = "jpeg", width = 12, height = 6, dpi = 300)

# ---- Class I heatmap --------------------------------------------------------------
aa_classI <- aa_all_top %>%
  filter(class == "ClassI", significant)

top_hm_t1 <- ggplot(aa_classI, aes(x = cluster, y = variant, fill = log10Pval)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "-log10(adj p)") +
  facet_wrap(~locus, scales = "free_y") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 10),
    panel.spacing = unit(0.3, "lines")
  ) +
  labs(
    title = "Per-Variant top HLA type I across Clusters",
    x = "Cluster",
    y = "Variant (position + residue)"
  )

ggsave(out_top_hm_t1, plot = top_hm_t1, device = "jpeg", width = 12, height = 10, dpi = 300)

# ----- Class II heatmap ------------------------------------------------------------------
aa_classII <- aa_all_top %>%
  filter(class == "ClassII", significant)

top_hm_t2 <- ggplot(aa_classII, aes(x = cluster, y = variant, fill = log10Pval)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "-log10(adj p)") +
  facet_wrap(~locus, scales = "free_y") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 10),
    panel.spacing = unit(0.3, "lines")
  ) +
  labs(
    title = "Per-Variant top HLA type II across Clusters",
    x = "Cluster",
    y = "Variant (position + residue)"
  )

ggsave(out_top_hm_t2, plot = top_hm_t2, device = "jpeg", width = 12, height = 12, dpi = 300)

# Save summary table
write.table(x = aa_all %>% select(-significant_cluster), file = output_summary, quote = F, sep="\t", row.names=F)