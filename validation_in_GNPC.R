# This script runs the replication in the GNPC cohort

library(tidyr)
library(purrr)
library(data.table)
library(stringr)
library(stringi)
library(fastmap)
library(readr)
library(RColorBrewer)
library(ggplot2)
library(WGCNA)
library(ggrepel)
library(dplyr)


# Loading module assignments for visit 1 from wgcna network object
bwnet_visit1_ref <- readRDS("data/visit1_results/bwnet_visit_1.rds")
bwnet_visit1_modassignments <- bwnet_visit1_ref$colors

# Loading expression data
# Loading pre-processed expression data(log2 transformed, batch corrected) for each visit where rows are proteins and columns are patients
counts_data_visit1 <- readRDS("data/visit1_results/counts_data_visit_1.rds")
counts_data_visit1 <- t(counts_data_visit1)

validation_data <- readRDS("data/gnpc_PD.rds")  # for cases
colnames(validation_data) <- gsub("_", ".", colnames(validation_data))

common_prots <-  Reduce(intersect, list(colnames(validation_data),
                         colnames(counts_data_visit1),
                         names(bwnet_visit1_modassignments)))

counts_data_visit1 <- counts_data_visit1[, common_prots]
validation_data <- validation_data[, common_prots, with = FALSE]
setcolorder(as.data.table(counts_data_visit1), common_prots)
setcolorder(validation_data, common_prots)
stopifnot(identical(colnames(validation_data), colnames(counts_data_visit1)))


# reindex module colors to EXACTLY match the column order used in the data
bwnet_visit1_modassignments <- bwnet_visit1_modassignments[common_prots]
if (anyNA(bwnet_visit1_modassignments)) {
  stop("Missing module labels for: ",
       paste(names(bwnet_visit1_modassignments)[is.na(bwnet_visit1_modassignments)], collapse = ", "))
}

stopifnot(identical(names(bwnet_visit1_modassignments), colnames(counts_data_visit1)))
stopifnot(identical(colnames(validation_data), colnames(counts_data_visit1)))

# Prepare multiData list (expression data)
multiData <- list(
  visit1 = list(data = counts_data_visit1),
  test = list(data = validation_data))

# Prepare multiColor list (module assignments for visit 1)
multiColor <- list(
  visit1 = bwnet_visit1_modassignments,
  test = bwnet_visit1_modassignments) # the module assignments for 'test' are unknown, 
                                      # but we still must provide a multiColor list of the same
                                      # length as multiData for modulePreservation() to run.
                                      
# NOTE: In modulePreservation(), the module color labels are primarily used for the
# reference network to define which modules to test. For the test/validation dataset,
# colors are not treated as "real" module assignments in the same way—it's common to
# supply a placeholder vector just to satisfy the required input.

# Run modulePreservation
mp_results <- modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  dataIsExpr = TRUE,
  referenceNetworks = 1,       
  testNetworks = 2, 
  nPermutations = 1000, 
  randomSeed = 1234,
  verbose = 3)
  
# Extract preservation stats
summary_stats_v1_test <- mp_results$preservation$Z[[1]][[2]]

# Convert to tidy data frame
plot_data_v1_test <- data.frame(
  moduleSize = summary_stats_v1_test$moduleSize,
  Zsummary   = summary_stats_v1_test$Zsummary.pres,
  comparison = "Visit1 (ref) vs Test",
  module     = rownames(summary_stats_v1_test)
)

plot_data <- plot_data_v1_test

saveRDS(mp_results, file = "data/mp_results.rds")
saveRDS(plot_data,  file = "data/plot_data.rds")


png("data/mpres_visit1_ref_vs_validation.png",
    width = 6000, height = 3600, res = 300)  

ggplot(plot_data, aes(x = moduleSize, y = Zsummary, color = comparison, label = module)) +
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
  


# create the median rank plot
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
 
summary_stats_v1_test <- mp_results$preservation$Z[[1]][[2]]
medianRank_v1_v2 <- get_medianRank(summary_stats_v1_test, mp_results, ref_idx = 1, test_idx = 2)

plot_data$medianRank <- medianRank_v1_v2
# filter gold module if needed
plot_data_filtered <- plot_data %>% dplyr::filter(module != "gold")


# --- Plot 2: medianRank vs moduleSize (lower = better) ---
png("mpres_medianRank_plot.png",
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

