# Creating consensus proteins clusters

library(readxl) 
library(ggplot2) 
library(WGCNA) 
library(gridExtra) 
library(clusterProfiler)
library(dplyr)

# Load files
# counts_data files for each visit should contain pre-processed(log2 transformed, batch corrected)
#  expression data for WGCNA where rows are proteins and columns are patients
counts_data_visit1 <- "data/visit1_results/counts_data_visit_1.rds" 
counts_data_visit2 <- "data/visit2_results/counts_data_visit_2.rds" 
counts_data_visit3 <- "data/visit3_results/counts_data_visit_3.rds" 
annotation_file <- "data/Somalogic_7-5K_protein_list.xlsx" # should contain annotations of all proteins
# (eg columns like SeqId, EntrezGeneID)

# Loading counts data for each visit
counts_data_visit1 <- readRDS(counts_data_visit1)
counts_data_visit1 <- t(counts_data_visit1)
counts_data_visit1[1:5,1:5]
dim(counts_data_visit1)

counts_data_visit2 <- readRDS(counts_data_visit2)
counts_data_visit2 <- t(counts_data_visit2)
counts_data_visit2[1:5,1:5]

counts_data_visit3 <- readRDS(counts_data_visit3)
counts_data_visit3 <- t(counts_data_visit3)
counts_data_visit3[1:5,1:5]

# Prepare multiData list (expression data for each visit)
multiData <- list(
  visit1 = list(data = counts_data_visit1),
  visit2 = list(data = counts_data_visit2),
  visit3 = list(data = counts_data_visit3))

# Check data for goodness
checkSets(multiData)

# Choose soft threshold power (do for all sets)
powers = c(1:10, seq(from = 12, to = 30, by = 2))
powerTables = vector(mode = "list", length = 3)

for (set in 1:3) {
  powerTables[[set]] = pickSoftThreshold(multiData[[set]]$data, 
                                        powerVector = powers,
                                        verbose = 2)
}

# Create a combined dataframe for all sets
sft.data.combined <- data.frame()
for (set in 1:3) {
  temp.data <- powerTables[[set]]$fitIndices
  temp.data$Dataset <- paste("Dataset", set)
  sft.data.combined <- rbind(sft.data.combined, temp.data)
}

# Visualization 1: Scale-free topology fit for all datasets
a1 <- ggplot(sft.data.combined, aes(Power, SFT.R.sq, label = Power, color = Dataset)) +
  geom_point() +
  geom_text(nudge_y = 0.02, size = 3) +
  geom_hline(yintercept = 0.8, color = 'red', linetype = "dashed") +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic() +
  facet_wrap(~Dataset, ncol = 3) +
  theme(legend.position = "none")

# Visualization 2: Mean connectivity for all datasets
a2 <- ggplot(sft.data.combined, aes(Power, mean.k., label = Power, color = Dataset)) +
  geom_point() +
  geom_text(nudge_y = 0.5, size = 3) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic() +
  facet_wrap(~Dataset, ncol = 3) +
  theme(legend.position = "none")

# Save the combined visualization
png("Multi_dataset_power_selection.png", width = 12, height = 8, units = "in", res = 300)
grid.arrange(a1, a2, nrow = 2)
dev.off()

# Creating consensus proteins clusters 
temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseConsensusModules(multiData,
                 maxBlockSize = 9000, 
                 corType = "bicor",
                 networkType = "signed",  # signed network better for biological interpretation
                 TOMType = "signed",  # TOM = Topological Overlap Measure
                 power = 18,   # choosing a power which has greater than 0.8 R sq value (above the red line) 
                               #  and has minimum mean connectivity
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE, # to have module eigen gene labels to be names of colours instead of numbers
                 randomSeed = 1234,
                 verbose = 3,
                 maxPOutliers = 0.07,
                 reassignThreshold = 0,    # avoid reassigning weak members
                 pamRespectsDendro = TRUE
                 )

cor <- temp_cor

module_eigengenes <- bwnet$MEs
table(bwnet$colors)

# Checking KEGG pathways for the clusters
all_pathway_results <- analyze_all_modules(
  bwnet = bwnet,  # WGCNA network object
  annotation_file = annotation_file,
  output_prefix = "consensus_cluster"
)
