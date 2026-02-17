library(tidyr)
library(purrr)
library(data.table)
library(writexl)
library(readxl)
library(stringr)
library(ggplot2)
library(WGCNA)
library(tibble)
library(KEGG.db)
library(gridExtra) 
library(clusterProfiler)
library(dplyr)

# Loading module assignments for each visit (as obtained from WGCNA on each visit. For visit 2 and visit 3,
# use renamed module assignments)
bwnet_visit1 <- readRDS("data/visit1_results/bwnet_visit_1.rds")
bwnet_visit1_modassignments <- bwnet_visit1$colors

visit2_modassignments_matched <- readRDS("data/visit2_results/visit2_modassignments_matched.rds") 

visit3_modassignments_matched_ref_visit2_modassignments_matched <- 
  readRDS("data/visit3_results/visit3_modassignments_matched_ref_visit2_modassignments_matched.rds") 

# Loading pre-processed expression data(log2 transformed, batch corrected) for each visit where rows are
#  proteins and columns are patients
counts_data_visit1 <- readRDS("data/visit1_results/counts_data_visit_1.rds")
counts_data_visit1 <- t(counts_data_visit1)
counts_data_visit1[1:5,1:5]

counts_data_visit2 <- readRDS("data/visit2_results/counts_data_visit_2.rds")
counts_data_visit2 <- t(counts_data_visit2)
counts_data_visit2[1:5,1:5]

counts_data_visit3 <- readRDS("data/visit3_results/counts_data_visit_3.rds")
counts_data_visit3 <- t(counts_data_visit3)
counts_data_visit3[1:5,1:5]

# Prepare multiData list (expression data for each visit)
multiData <- list(
  visit1 = list(data = counts_data_visit1),
  visit2 = list(data = counts_data_visit2),
  visit3 = list(data = counts_data_visit3))

# Prepare multiColor list (module assignments for each visit)
multiColor <- list(
  visit1 = bwnet_visit1_modassignments,
  visit2 = visit2_modassignments_matched, 
  visit3 = visit3_modassignments_matched_ref_visit2_modassignments_matched)

# Run modulePreservation - both visit 1 and 2 are given as reference
mp_results <- modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  dataIsExpr = TRUE,
  maxGoldModuleSize = 1000,
  maxModuleSize = 6000, 
  referenceNetworks = c(1, 2),        # Ref: visit1 and visit2
  testNetworks = list(                # Must align with each reference
    c(2, 3),                          # Tests for visit1 reference: visits 2 & 3
    c(3)                           # Tests for visit2 reference: visit 3
  ), 
  nPermutations = 1000,
  randomSeed = 1234,
  verbose = 3)

# Extract preservation stats
summary_stats_v1_v2 <- mp_results$preservation$Z[[1]][[2]]
summary_stats_v1_v3 <- mp_results$preservation$Z[[1]][[3]]
summary_stats_v2_v3 <- mp_results$preservation$Z[[2]][[3]]

# Convert to tidy data frames
plot_data_v1_v2 <- data.frame(
  moduleSize = summary_stats_v1_v2$moduleSize,
  Zsummary   = summary_stats_v1_v2$Zsummary.pres,
  comparison = "Visit1 (ref) vs Visit2",
  module     = rownames(summary_stats_v1_v2)
)

plot_data_v1_v3 <- data.frame(
  moduleSize = summary_stats_v1_v3$moduleSize,
  Zsummary   = summary_stats_v1_v3$Zsummary.pres,
  comparison = "Visit1 (ref) vs Visit3",
  module     = rownames(summary_stats_v1_v3)
)

plot_data_v2_v3 <- data.frame(
  moduleSize = summary_stats_v2_v3$moduleSize,
  Zsummary   = summary_stats_v2_v3$Zsummary.pres,
  comparison = "Visit2 (ref) vs Visit3",
  module     = rownames(summary_stats_v2_v3)
)

# Bind all together
plot_data <- bind_rows(
  plot_data_v1_v2,
  plot_data_v1_v3,
  plot_data_v2_v3
)

saveRDS(mp_results, file = "data/mp_results_module_pres.rds")
saveRDS(plot_data, file = "data/plot_data_module_pres.rds")

png("data/module_pres_plot.png",
    width = 6000, height = 3600, res = 300)  

# Remove gold module before plotting
filtered_plot_data <- plot_data %>%
filter(module != 'gold')

ggplot(filtered_plot_data, aes(x = moduleSize, y = Zsummary, color = comparison, label = module)) +
  geom_point(size = 2.4, alpha = 0.85) +   
  geom_text_repel(size = 5.4, box.padding = 0.5, point.padding = 0.2,
                  show.legend = FALSE, max.overlaps = 200) +
  geom_hline(yintercept = c(2, 10), linetype = "dashed", color = "red") +

  # use original (linear) scales with padding + pretty ticks
  scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                     labels = scales::label_number(big.mark = ","),
                     expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 8),
                     expand = expansion(mult = c(0.02, 0.06))) +

  labs(x = "Module Size", y = "Zsummary") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 16),
    axis.text   = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title= element_text(size = 15)
  ) +
  coord_cartesian(clip = "off")  # lets repelled labels hang outside the panel if needed

dev.off()

# ##### Extracting medianRank and plotting it 
# # --- extract medianRank robustly and add to plot_data ---

# helper to safely pull medianRank.pres from a stats object or from observed
get_medianRank <- function(summary_obj, mp_results, ref_idx, test_idx) {
  # try direct column in the summary (what you used for Zsummary)
  if (!is.null(summary_obj$medianRank.pres)) {
    return(summary_obj$medianRank.pres)
  }
  # fallback: try observed slot where medianRank may live
  # mp_results$preservation$observed is indexed similarly to Z: [[ref]][[test]]
  observed_slot <- tryCatch(mp_results$preservation$observed[[ref_idx]][[test_idx]],
                            error = function(e) NULL)
  if (!is.null(observed_slot) && !is.null(observed_slot$medianRank.pres)) {
    return(observed_slot$medianRank.pres)
  }
  # if nothing found, return NA vector (same length as moduleSize if available)
  len <- if (!is.null(summary_obj$moduleSize)) length(summary_obj$moduleSize) else 0
  return(rep(NA_real_, len))
}

# extract medianRank vectors
medianRank_v1_v2 <- get_medianRank(summary_stats_v1_v2, mp_results, ref_idx = 1, test_idx = 2)
medianRank_v1_v3 <- get_medianRank(summary_stats_v1_v3, mp_results, ref_idx = 1, test_idx = 3)
medianRank_v2_v3 <- get_medianRank(summary_stats_v2_v3, mp_results, ref_idx = 2, test_idx = 3)

# add to the individual data.frames
plot_data_v1_v2$medianRank <- medianRank_v1_v2
plot_data_v1_v3$medianRank <- medianRank_v1_v3
plot_data_v2_v3$medianRank <- medianRank_v2_v3

# bind again
plot_data <- bind_rows(
  plot_data_v1_v2,
  plot_data_v1_v3,
  plot_data_v2_v3
)

View(plot_data)

# remove 'gold' module 
plot_data_filtered <- plot_data %>% filter(module != "gold")
View(plot_data_filtered)

# --- Plot : medianRank vs moduleSize (lower = better) ---

png("data/mpres_medianRank_plot_long.png",
    width = 3000, height = 4800, res = 300) 
  ggplot(plot_data_filtered, aes(x = moduleSize, y = medianRank, color = comparison, label = module)) +
    geom_point(size = 2.4, alpha = 0.85) +
    geom_text_repel(size = 5.4, box.padding = 0.5, point.padding = 0.2,
                    show.legend = FALSE, max.overlaps = 200) +
    # smaller medianRank is better -> reverse the y-axis for intuitive plotting
    scale_y_reverse(breaks = scales::breaks_pretty(n = 8),
                    expand = expansion(mult = c(0.02, 0.06))) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                        labels = scales::label_number(big.mark = ","),
                        expand = expansion(mult = c(0.02, 0.06))) +
    labs(x = "Module Size", y = "medianRank") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom")
dev.off()

