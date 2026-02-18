# This script performs data preprocessing and quality control

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
library(sva)
library(GOSemSim)
library(org.Hs.eg.db)
library(EWCE)
library(grid)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(ggalluvial)

library(dplyr)
library(KEGG.db)

# Load custom functions
source("all_custom_functions.R")


# Loading expression and phenotype data
my_adat_visit <- read_adat("Tracking_proteomics.adat") # expression data
phenotype_data <- readRDS("Tracking_prot_clin.Rds") # phenotype data

# Subsetting data for 3 visits
### Note - The visit_num and corresponding visits in the data are 1(visit1),4(visit2),7(visit3),9(visit4),10(visit5),11(visit6)

# converts the my_adat_visit object to a data frame and then filters the rows of the data frame
# to include only those rows where the SampleId column ends with _v1, _v4, or _v7
my_adat_visits_filtered <- as.data.frame(my_adat_visit) %>%
filter(grepl("_v[147]$", SampleId))
dim(my_adat_visits_filtered) 

# Remove duplicates from my_adat_visits_filtered
my_adat_visits_filtered_rm_dup <- my_adat_visits_filtered[!duplicated(my_adat_visits_filtered$SampleId), ]
dim(my_adat_visits_filtered_rm_dup) 

# removing controls and keeping passed samples only
cleanData <- my_adat_visits_filtered_rm_dup %>%
filter(SampleType == "Sample") %>%                # rm control samples
filter(RowCheck == "PASS")                        # # keep passed samples only

# removing undesired columns
filtered_cleanData <- cleanData[, !names(cleanData) %in% c( "PlateId","PlateRunDate","ScannerID","PlatePosition", 
"SlideId", "Subarray", "SampleType", "PercentDilution", "SampleMatrix", "Barcode", "Barcode2d", "SampleName", "SampleNotes", 
"AliquotingNotes", "SampleDescription", "AssayNotes", "TimePoint", "ExtIdentifier", "SsfExtId", "SampleGroup", "SiteId", 
"SubjectID", "CLI", "RMA", "HybControlNormScale", "RowCheck", "NormScale_20", "NormScale_0_005", "NormScale_0_5", "ANMLFractionUsed_20",
 "ANMLFractionUsed_0_005", "ANMLFractionUsed_20", "ANMLFractionUsed_0_005", "ANMLFractionUsed_0_5" )]

sum(is.na(filtered_cleanData)) # check missing data

# Combine the expression data with phenotype data in one df
Data <- left_join(filtered_cleanData, phenotype_data,  by = "SampleId")

# Subsetting the data to include common sample names across the first three visits (visit numbers 1, 4, and 7),
# such that the sample size for each of the three visits is equal
# Extract visit numbers from SampleId
visits <- unique(str_extract(Data$SampleId, "_v[1-7]+$"))

## Filter for visits 1, 4 and 7
# Count samples for each visit
visit_counts <- Data %>%
  mutate(SampleName = str_remove(SampleId, "_v[1-7]+$")) %>%
  filter(str_extract(SampleId, "_v[1-7]+$") %in% visits) %>%
  count(SampleName, name = "n") %>%
  filter(n == 3)

# Filter dataframe for common samples across visits 1, 4, 7
common_samples_df <- Data %>%
  filter(str_remove(SampleId, "_v[1-7]+$") %in% visit_counts$SampleName)

# Check sample size for each visit
common_samples_df %>%
  mutate(Visit = str_extract(SampleId, "_v[1-7]+$")) %>%
  count(Visit)

# we get this sample size for each visit
#   Visit   n
# 1   _v1 794
# 2   _v4 794
# 3   _v7 794
saveRDS(common_samples_df, file = "common_samples_df.rds")
dim(common_samples_df) 



## Filtering out non human proteins from common_samples_df
Annotation_data <- read_excel("Somalogic_7-5K_protein_list.xlsx")
View(Annotation_data)

# Modifying the seqID column
Annotation_data <- Annotation_data %>%
  mutate(SeqId = gsub("^(\\d+)-(\\d+)$", "seq.\\1.\\2", SeqId))

# Modifying the dataframe to have just SeqId and EntrezGeneID
Annotation_data_new <- Annotation_data %>%
filter(Organism == "Human") %>% 
filter(Type == "Protein")

dim(Annotation_data_new) 
View(Annotation_data_new)

# Extract the SeqId values
seq_ids <- Annotation_data_new$SeqId

# Keep only columns in common_samples_df that are also in seq_ids
common_samples_df <- common_samples_df %>% 
  select(any_of(seq_ids))


# PCA of raw expression data and coloring with technical variables (eg; Site, Aliquot, Wells and Visit) to check for batch effects
common_samples_df_expr <- common_samples_df[, is_seq(names(common_samples_df))] # selecting expr data
PCA_common_samples_df <- common_samples_df_expr
dim(PCA_common_samples_df) #2382 7288
sum(is.na(PCA_common_samples_df))

my_pca <- prcomp(PCA_common_samples_df)

pca.dat <- my_pca$x

pca.var <- my_pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)


png("PCA_common_samples_df.png")
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %')) +
  ggtitle("PCA 3 visit log2 data")
dev.off()

png("PCA 3 visit batch effect check.png")
par(mfrow = c(2, 2))

title = "PCA plot 1 - Coloured wrt common_samples_df$Site.y"  

# Create a scatter plot of the first and second principal components
plot(pca.dat[,1],
     pca.dat[,2],
     pch=15,          
     bty='l',        
     xlab='PCA 1',    
     ylab='PCA 2',  
     cex=0.8,         
     cex.axis=1,     
     cex.main=1,      
     cex.lab=1.2,    
     main = title,    
     xlim=c(min(pca.dat[,1]-abs(0.2*min(pca.dat[,1]))),      
            max(pca.dat[,1]+abs(0.2*max(pca.dat[,1])))),
     ylim=c(min(pca.dat[,2]-abs(0.2*min(pca.dat[,2]))),
            max(pca.dat[,2]+abs(0.2*max(pca.dat[,2])))) ,     
     col= brewer.pal(name = "Set3", n = 12)[as.factor(common_samples_df$Site.y)])  



title = "PCA plot 2 - Coloured wrt common_samples_df$Aliquot"   

# Create a scatter plot of the first and second principal components
plot(pca.dat[,1],
     pca.dat[,2],
     pch=15,         
     bty='l',         
     xlab='PCA 1',    
     ylab='PCA 2',   
     cex=0.8,         
     cex.axis=1,      
     cex.main=1,      
     cex.lab=1.2,     
     main = title,    
     xlim=c(min(pca.dat[,1]-abs(0.2*min(pca.dat[,1]))),      
            max(pca.dat[,1]+abs(0.2*max(pca.dat[,1])))),
     ylim=c(min(pca.dat[,2]-abs(0.2*min(pca.dat[,2]))),
            max(pca.dat[,2]+abs(0.2*max(pca.dat[,2])))) ,
     col= brewer.pal(name = "Set3", n = 12)[as.factor(common_samples_df$Aliquot)]) 


title = "PCA plot 3 - Coloured wrt common_samples_df$Visit" 

# Create a scatter plot of the first and second principal components
plot(pca.dat[,1],
     pca.dat[,2],
     pch=15,          
     bty='l',         
     xlab='PCA 1',    
     ylab='PCA 2',    
     cex=0.8,         
     cex.axis=1,      
     cex.main=1,      
     cex.lab=1.2,     
     main = title,    
     xlim=c(min(pca.dat[,1]-abs(0.2*min(pca.dat[,1]))),      
            max(pca.dat[,1]+abs(0.2*max(pca.dat[,1])))),
     ylim=c(min(pca.dat[,2]-abs(0.2*min(pca.dat[,2]))),
            max(pca.dat[,2]+abs(0.2*max(pca.dat[,2])))) ,
     col= brewer.pal(name = "Set3", n = 12)[as.factor(common_samples_df$Visit)])  



title = "PCA plot 4 - Coloured wrt common_samples_df$Wells"   

# Create a scatter plot of the first and second principal components
plot(pca.dat[,1],
     pca.dat[,2],
     pch=15,          
     bty='l',         
     xlab='PCA 1',    
     ylab='PCA 2',   
     cex=0.8,         
     cex.axis=1,      
     cex.main=1,      
     cex.lab=1.2,     
     main = title,    
     xlim=c(min(pca.dat[,1]-abs(0.2*min(pca.dat[,1]))),     
            max(pca.dat[,1]+abs(0.2*max(pca.dat[,1])))),
     ylim=c(min(pca.dat[,2]-abs(0.2*min(pca.dat[,2]))),
            max(pca.dat[,2]+abs(0.2*max(pca.dat[,2])))) ,
     col= brewer.pal(name = "Set3", n = 12)[as.factor(common_samples_df$Wells)])

dev.off()

## PCA revealed batch effects due to Site 
# Using sva package for batch correction
# log2 transformation of expression data 
common_samples_df_expr <- common_samples_df_expr %>% log2()

# Batch correction
batch = common_samples_df$Site.y  # initializing batch as Site.y as we need to correct data for sites
sum(is.na(common_samples_df$Site.y )) 

param <- MulticoreParam(workers = 8)

Combat_corrected_common_samples_df_expr =  ComBat(dat= t(common_samples_df_expr), batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE, BPPARAM = param)

sum(is.na(Combat_corrected_common_samples_df_expr))

Combat_corrected_common_samples_df_expr <- t(Combat_corrected_common_samples_df_expr)

## Note: Run PCA again on batch corrected expression data to check if the batch effect is removed


# combine log transformed and batch corrected expression data with phenotype data in one df
stopifnot(nrow(common_samples_df_expr) == nrow(phenotype_data))

all_data <- bind_cols(
    common_samples_df_expr,
    phenotype_data)

# drop rows which have missing age_at_visit and gender values
all_data <- all_data[!is.na(all_data$age_at_visit), ] # all missing gender rows are eliminated as well here

saveRDS(all_data, file = "data/preprocessed_expr_n_pheno_data_3_visits.rds")