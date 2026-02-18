library(tidyr)
library(purrr)
library(data.table)
library(writexl)
library(readxl)
library(stringr)
library(stringi)
library(fastmap)
library(readr)
library(RColorBrewer)
library(SomaDataIO)
library(ggplot2)
library(WGCNA)
library(CorLevelPlot)
library(gridExtra) 
library(tibble)
library(writexl)
library(sva)
library(clusterProfiler)
library(GOSemSim)
library(org.Hs.eg.db)
library(EWCE)
library(grid)
library(enrichplot)
library(ggrepel)
library(ggalluvial)

library(dplyr)
library(KEGG.db)

# Load custom functions
source("all_custom_functions.R")

#### 1. Preparing data for WGCNA workflow -----------------------------------------------------------
# Load preprocessed data for 3 visits, as obtained from the preprocessing script
all_data <- readRDS("data/preprocessed_expr_n_pheno_data_3_visits.rds")

# Partitioning data by visit 
visit_1_df <- all_data %>% 
filter(grepl("_v1$",SampleId))
dim(visit_1_df) 

visit_2_df <- all_data %>% 
filter(grepl("_v4",SampleId))
dim(visit_2_df) 

visit_3_df <- all_data %>% 
filter(grepl("_v7$",SampleId))
dim(visit_3_df) 

# Sample outliers were assessed for each visit independently using two methods: hierarchical clustering and principal component analysis
# repeat these methods for each visit
visit_1_df_expr <- visit_1_df[, is_seq(names(visit_1_df))]  # expression data
htree <- hclust(dist(visit_1_df_expr), method = "average") # input data rows should be samples and columns should be genes/proteins
png("data/visit1_results/hclust_plot.png", width = 3500, height = 2000, res = 300)
plot(htree) 
dev.off()

pca <- prcomp(visit_1_df_expr)  # rows should be samples and columns should be genes/proteins
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

png("data/visit3_results/pca_plot.png", width = 2500, height = 3000, res = 300)
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
dev.off()


#### 2. WGCNA workflow (entire workflow to be repeated for each visit data)-----------------------------------------------------------
# results saved from this code block will be used for other analyses in this script and in other scripts.
# uncomment the code below based on the visit number you're running wgcna for
# For visit 1
visit_1_df <- visit_1_df
#For visit 2 
# visit_1_df <- visit_2_df
#For visit 3 
# visit_1_df <- visit_3_df

View(visit_1_df)
visit_1_df[1:5, 1:5]

# remove outlier samples as identified above, in the script. Any sample found as outlier in one visit is removed from all other visits as
# well, to keep the samples and sample size consistent across visits. 
visit_1_df <- visit_1_df %>%
  filter(!SampleId %in% samples_to_remove) # samples_to_remove variable should contain the outlier samples for the visit wgcna is ran for
dim(visit_1_df)

## Preparing counts_data (expression data in the format for wgcna)
visit_1_df_expr <- visit_1_df[, is_seq(names(visit_1_df))]
dim(visit_1_df_expr) 
sum(is.na(visit_1_df_expr)) 
sum(visit_1_df_expr < 0) 

# Creating rownames (which are samples) from the SampleId column
rownames(visit_1_df_expr) <- visit_1_df$SampleId
# Transposing dataframe to have protein names as rownames and sample names as columns names
counts_data <- t(visit_1_df_expr)
counts_data[1:5,1:5]
dim(counts_data)

saveRDS(counts_data, file = "data/visit1_results/counts_data_visit_1.rds")

## Preparing phenodata (phenotype data in the format for wgcna)
pheno_data <- visit_1_df %>%
  select("moca_total", "moca_total_adj" , "seman_flu_score", "updrs_i", "updrs_ii", "updrs_iii", "updrs_iv","disease_duration_diag", "age_at_visit", "pdq8_total", "ess_total","eq5d_index",
 "rbd_total", "gcsi_total", "leeds_anx_total", "hoehn_and_yahr_stage",
"nmss_sleep","ldopa_bin","orthostatic_hypotension","quip_sex" , "quip_med")
rownames(pheno_data) <- visit_1_df$SampleId
dim(pheno_data)

saveRDS(pheno_data, file = "data/visit1_results/pheno_data_visit_1.rds")


all(rownames(pheno_data) %in% colnames(counts_data)) 
all(rownames(pheno_data) == colnames(counts_data))

# Outlier detection
gsg <- goodSamplesGenes(t(counts_data))  # rows should be samples and columns should be genes
summary(gsg)
gsg$allOK 

table(gsg$goodGenes) 
table(gsg$goodSamples) 

# Preparing wgcna format specific input expr data
counts_input_data <- as.data.frame(t(counts_data))
counts_input_data[1:5,1:5]
dim(counts_input_data) 

# Network Construction 
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function - function in WGCNA package for picking soft threshold
# Choose the soft threshold that produces high similarity and scale free network
sft <- pickSoftThreshold(counts_input_data,, #
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)

sft.data <- sft$fitIndices

# Visualization to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
  

# Ideally choose a power which has greater than 0.8 R sq value (above the red line) and has minimum mean connectivity
# If more than one powers are above the red line, choose the one which has greater than 0.8 R sq value but not excessively large 
png("data/visit1_results/Picking_soft_threshold_power.png", width = 2000, height = 2500, res = 300)
grid.arrange(a1, a2, nrow = 2)
dev.off()

# Convert matrix to numeric
counts_input_data[] <- sapply(counts_input_data, as.numeric)

soft_power <- 18  # we found 18 to be a suitable power
temp_cor <- cor
cor <- WGCNA::cor

# Creating adjacency matrix, weighted correlation matrix, calculating proximity measures(that is topological overlap matrix), creating hierarchical clustering, merging modules that are simialr or have a high correlation - all this done by one function blockwiseModules()
bwnet <- blockwiseModules(counts_input_data,
                 maxBlockSize = 9000,
                 corType = "bicor",
                 networkType = "signed", # signed network better for biological interpretation - eg a module and all genes within it has negative correlation with a trait,if the correlation sign is negative value 
                 TOMType = "signed",  # TOM = Topological Overlap Measure
                 power = soft_power,
                 mergeCutHeight = 0.25,#0.25, #threshold to merge simialr modules
                 numericLabels = FALSE, # to have module eigen gene labels to be names of colours instead of numbers
                 randomSeed = 1234,
                 verbose = 3,
                 maxPOutliers = 0.07,
                 reassignThreshold = 0,    # avoid reassigning weak members
                 pamRespectsDendro = TRUE
                 )


cor <- temp_cor

# Plot the dendrogram and the module colors before and after merging 
png("data/visit1_results/cluster_dendogram.png", width = 2000, height = 1500, res = 300)
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)
dev.off()

saveRDS(bwnet, file = "data/visit1_results/bwnet_visit_1.rds")

module_eigengenes <- bwnet$MEs
saveRDS(module_eigengenes, file = "data/visit1_results/module_eigengenes_visit_1.rds")

#### end of WGCNA workflow (repeat for next visit data)---


#### 3. Renaming the modules of visit2 and visit3 -------------------------------------------------------------------------
# loading module assignments from all visits as generated from the above steps
bwnet_visit1 <- readRDS("data/visit1_results/bwnet_visit_1.rds")
bwnet_visit1_modassignments <- bwnet_visit1$colors
head(bwnet_visit1_modassignments)
table(bwnet_visit1_modassignments)

bwnet_visit2 <- readRDS("data/visit2_results/bwnet_visit_2.rds")
bwnet_visit2_modassignments <- bwnet_visit2$colors
head(bwnet_visit2_modassignments)
table(bwnet_visit2_modassignments)

bwnet_visit3 <- readRDS("data/visit3_results/bwnet_visit_3.rds")
bwnet_visit3_modassignments <- bwnet_visit3$colors
head(bwnet_visit3_modassignments)
table(bwnet_visit3_modassignments)

# Renaming visit2 module assignments wrt visit 1
visit2_modassignments_matched <- matchLabels(bwnet_visit2_modassignments, bwnet_visit1_modassignments)
table(bwnet_visit1_modassignments,visit2_modassignments_matched)
table(bwnet_visit2_modassignments,visit2_modassignments_matched)
table(visit2_modassignments_matched)
saveRDS(visit2_modassignments_matched, file = "data/visit2_results/visit2_modassignments_matched.rds")

# Renaming visit3 modules assignments wrt to visit2_modassignments_matched 
visit3_modassignments_matched_ref_visit2_modassignments_matched <- matchLabels(bwnet_visit3_modassignments, visit2_modassignments_matched)
table(visit2_modassignments_matched, visit3_modassignments_matched_ref_visit2_modassignments_matched)
table(bwnet_visit1_modassignments,visit3_modassignments_matched_ref_visit2_modassignments_matched )
table(bwnet_visit3_modassignments,visit3_modassignments_matched_ref_visit2_modassignments_matched )
table(visit3_modassignments_matched_ref_visit2_modassignments_matched)
saveRDS(visit3_modassignments_matched_ref_visit2_modassignments_matched, file = "data/visit3_results/visit3_modassignments_matched_ref_visit2_modassignments_matched.rds")

# You can rename visit3 modules assignments wrt to bwnet_visit1_modassignments
# To see if the renamed module assignments of visit 3 change based on the reference being visit 2 or visit 1 

# visit3_modassignments_matched_ref_visit1 <- matchLabels(bwnet_visit3_modassignments, bwnet_visit1_modassignments)
# table(visit2_modassignments_matched, visit3_modassignments_matched_ref_visit1)
# table(bwnet_visit1_modassignments,visit3_modassignments_matched_ref_visit1)
# table(bwnet_visit3_modassignments,visit3_modassignments_matched_ref_visit1)
# table(visit3_modassignments_matched_ref_visit2_modassignments_matched,visit3_modassignments_matched_ref_visit1)
# table(visit3_modassignments_matched_ref_visit1)

#####  WE FOUND NO MATTER THE REFERENCE WE CHOOSE, WE GET THE SAME RENAMED MODULE ASSIGNMENTS ->

# > table(visit3_modassignments_matched_ref_visit2_modassignments_matched)
# visit3_modassignments_matched_ref_visit2_modassignments_matched
#     black      blue     brown     green      grey   magenta       red turquoise 
#        25       777       202        45       985        30        23      5133 
#    yellow 
#        68 
# > table(visit3_modassignments_matched_ref_visit1)
# visit3_modassignments_matched_ref_visit1
#     black      blue     brown     green      grey   magenta       red turquoise 
#        25       777       202        45       985        30        23      5133 
#    yellow 
#        68 

### ------



#### 4. Creating and saving combined supplementary table with RELABELED modules for all visits ----------------------------------------
# Load annotation data for all proteins
Annotation_data <- read_excel("Somalogic_7-5K_protein_list.xlsx") %>%
  select("SeqId", "EntrezGeneSymbol") %>%
  mutate(SeqId = gsub("^(\\d+)-(\\d+)$", "seq.\\1.\\2", SeqId))

# Load module assignment vectors
bwnet_visit1 <- readRDS("data/visit1_results/bwnet_visit_1.rds")
bwnet_visit1_modassignments <- bwnet_visit1$colors
names(bwnet_visit1_modassignments)[1:10]

visit2_modassignments_matched <- readRDS("data/visit2_results/visit2_modassignments_matched.rds")
names(visit2_modassignments_matched) <- names(bwnet_visit1_modassignments)

visit3_modassignments_matched_ref_visit2_modassignments_matched <- readRDS("data/visit3_results/visit3_modassignments_matched_ref_visit2_modassignments_matched.rds")
names(visit3_modassignments_matched_ref_visit2_modassignments_matched) <- names(bwnet_visit1_modassignments)

supp_table_all_visits <- data.frame(
  Seq_ID = names(bwnet_visit1$colors),
  Visit1_Module = bwnet_visit1$colors,  # Original (reference)
  Visit2_Module = visit2_modassignments_matched,  # Relabeled
  Visit3_Module = visit3_modassignments_matched_ref_visit2_modassignments_matched,  # Relabeled
  stringsAsFactors = FALSE
)
# Merge with supplementary table
supp_table_all_visits <- supp_table_all_visits %>%
  left_join(Annotation_data, by = c("Seq_ID" = "SeqId"))

# Reorder columns for better readability (optional)
supp_table_all_visits <- supp_table_all_visits %>%
  select(Seq_ID, EntrezGeneSymbol, EntrezGeneID, 
         Visit1_Module, Visit2_Module, Visit3_Module)

# View first few rows
head(supp_table_all_visits)
dim(supp_table_all_visits)

table(supp_table_all_visits$Visit1_Module)
table(supp_table_all_visits$Visit2_Module)
table(supp_table_all_visits$Visit3_Module)

head(table(supp_table_all_visits))

# removing grey module from the supplementary table
supp_table_common_grey_removed <- supp_table_all_visits %>%
  filter(!(Visit1_Module == "grey" & 
           Visit2_Module == "grey" & 
           Visit3_Module == "grey"))
dim(supp_table_common_grey_removed)
head(supp_table_common_grey_removed)

supp_table_common_grey_removedandNA <- supp_table_common_grey_removed %>%
  mutate(Visit1_Module = na_if(Visit1_Module, "grey"),
         Visit2_Module = na_if(Visit2_Module, "grey"),
         Visit3_Module = na_if(Visit3_Module, "grey"))

# Check module counts to verify grey is gone
table(supp_table_common_grey_removedandNA$Visit1_Module)
table(supp_table_common_grey_removedandNA$Visit2_Module)
table(supp_table_common_grey_removedandNA$Visit3_Module)

write_xlsx(supp_table_common_grey_removedandNA, 
           "data/Supplementary_Table_Module_Assignments_All_Visits_Annotated_AllGrey_removed.xlsx")

#### ---



#### 5. WGCNA Correlation plot between module eigengenes and clinical phenotypes (entire code chunk to be repeated for each visit data)---------------------------------------------------------------
## We would have saved module eigengenes for each visit but we need to create module eigengenes following the creation of 
# relabelled module assignments for visit 2 and visit 3. This ensures the module names are consistent between the visits

# Load module assignment vectors
bwnet_visit1 <- readRDS("data/visit1_results/bwnet_visit_1.rds")
bwnet_visit1_modassignments <- bwnet_visit1$colors
names(bwnet_visit1_modassignments)[1:10]

visit2_modassignments_matched <- readRDS("data/visit2_results/visit2_modassignments_matched.rds")
names(visit2_modassignments_matched) <- names(bwnet_visit1_modassignments)

visit3_modassignments_matched_ref_visit2_modassignments_matched <- readRDS("data/visit3_results/visit3_modassignments_matched_ref_visit2_modassignments_matched.rds")
names(visit3_modassignments_matched_ref_visit2_modassignments_matched) <- names(bwnet_visit1_modassignments)

# Prepare module eigengenes for each visit

MEs_v1 <- readRDS("data/visit1_results/module_eigengenes_visit_1.rds") # visit 1 modules eigengenes can be loaded as it was saved above

datExpr_v2 <- readRDS("data/visit2_results/counts_data_visit_2.rds")     # visit 2 expression matrix used for WGCNA 
MEs_v2     <- moduleEigengenes(datExpr_v2, visit2_modassignments_matched)$eigengenes

datExpr_v3 <- readRDS("data/visit3_results/counts_data_visit_3.rds")     # expression matrix used for WGCNA    # expression matrix used for WGCNA
MEs_v3     <- moduleEigengenes(datExpr_v3, visit3_modassignments_matched_ref_visit2_modassignments_matched)$eigengenes

# Load phenotype data for all visits (which were generated above - before doing wgcna)
pheno_data_visit_1 <- readRDS("data/visit1_results/pheno_data_visit_1.rds")
pheno_data_visit_2 <- readRDS("data/visit2_results/pheno_data_visit_2.rds")
pheno_data_visit_3 <- readRDS("data/visit3_results/pheno_data_visit_3.rds")

# visualize module-trait association as a heatmap (repeat the entire below chunk of code for each visit module eigengene and phenotype data)
heatmap.data <- merge(MEs_v1, pheno_data_visit_1, by = 'row.names') 
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

heatmap.data <- heatmap.data %>%
  select(-c(MEgrey)) 

# CHANGE the MODULE NAMES based on  DIFFERENT VISITS (some visits have more or less modules than others)
cols_to_keep <- c("MEbrown","MEblue","MEturquoise","MEgreen","MEyellow","MEblack", "MEred", "moca_total", "moca_total_adj" , "seman_flu_score", "updrs_i", "updrs_ii", "updrs_iii", "updrs_iv","disease_duration_diag", "age_at_visit", "pdq8_total", "ess_total","eq5d_index",
 "rbd_total", "gcsi_total", "leeds_anx_total", "hoehn_and_yahr_stage",
"nmss_sleep","ldopa_bin","orthostatic_hypotension","quip_sex" , "quip_med")

heatmap.data <- heatmap.data %>%
  select(all_of(cols_to_keep))

names(heatmap.data)

# creating correlation plot
png("data/visit1_results/correlation_plot_visit1.png",width = 4000, height = 1500, res = 300, bg = "white")
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[8:21], # should include phenotypes
             y = names(heatmap.data)[1:7], # should include module names 
             col = c("blue1", "skyblue", "white", "pink", "red"),
              rotLabX = 45)
dev.off()

## ---