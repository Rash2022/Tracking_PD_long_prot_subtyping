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
library(GOSemSim)
library(org.Hs.eg.db)
library(EWCE)
library(grid)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(ggalluvial)
library(enrichR)

library(dplyr)
library(KEGG.db)

# Load custom functions
source("all_custom_functions_code_git.R")

# There are 8 independent analyses in the script. You can run each analysis separately.

#### 1. Quick Pathway enrichment check (for all modules of one visit)------------------------------------------------------------------------
# # Repeat with different input data corresponding to different visits
all_pathway_results <- analyze_all_modules(
  bwnet = bwnet,  # Your WGCNA network object of a particular visit
  annotation_file = "Somalogic_7-5K_protein_list.xlsx",
  output_prefix = "visit3" # change for different visits 
)
### ----

#### 2. Publication ready Pathway enrichment facet plot for modules across visits ------------------------------------------------------------

# ============================================================================
# STEP 1: Load annotation data and module assignments
# ============================================================================

Annotation_data <- read_excel("data/Somalogic_7-5K_protein_list.xlsx") %>%
  select("SeqId", "EntrezGeneID", "EntrezGeneSymbol") %>%
  mutate(SeqId = gsub("^(\\d+)-(\\d+)$", "seq.\\1.\\2", SeqId))

# Load module assignment vectors
bwnet_visit1 <- readRDS("data/visit1_results/bwnet_visit_1.rds")
bwnet_visit1_modassignments <- bwnet_visit1$colors
names(bwnet_visit1_modassignments)[1:10]

visit2_modassignments_matched <- readRDS("data/visit2_results/visit2_modassignments_matched.rds")
names(visit2_modassignments_matched) <- names(bwnet_visit1_modassignments)

visit3_modassignments_matched_ref_visit2_modassignments_matched <- readRDS("data/visit3_results/visit3_modassignments_matched_ref_visit2_modassignments_matched.rds")
names(visit3_modassignments_matched_ref_visit2_modassignments_matched) <- names(bwnet_visit1_modassignments)

# ============================================================================
# STEP 2: Function to extract genes from modules (takes module assignments directly)
# ============================================================================

extract_module_genes <- function(module_assignments, module_color, annotation_data) {
  # Ensure module_assignments is a named vector
  if (is.null(names(module_assignments))) {
    stop("Error: 'module_assignments' must be a named vector with SeqIds as names.")
  }
  
  # Get SeqIds belonging to this module
  module_proteins <- names(module_assignments)[module_assignments == module_color]
  
  # Convert to dataframe for joining
  module_proteins_df <- data.frame(SeqId = module_proteins, stringsAsFactors = FALSE)
  
  # Map to annotation data
  associated_proteins <- left_join(module_proteins_df, annotation_data, by = "SeqId")
  
  # Return valid EntrezGene IDs
  entrez_genes <- na.omit(associated_proteins$EntrezGeneID)
  
  return(entrez_genes)
}

# ============================================================================
# STEP 3: Extract genes for all modules and visits
# ============================================================================

modules_of_interest <- c("black","blue","brown","turquoise","red")

gene_lists <- list()

visit1_modules <- bwnet_visit1_modassignments
visit2_modules <- visit2_modassignments_matched
visit3_modules <- visit3_modassignments_matched_ref_visit2_modassignments_matched

for (module in modules_of_interest) {
  
  # Visit 1
  if (module %in% unique(visit1_modules)) {
    genes_v1 <- extract_module_genes(visit1_modules, module, Annotation_data)
    if (length(genes_v1) > 0) {
      gene_lists[[paste0("Visit1_", str_to_title(module))]] <- genes_v1
    }
  }
  
  # Visit 2
  if (module %in% unique(visit2_modules)) {
    genes_v2 <- extract_module_genes(visit2_modules, module, Annotation_data)
    if (length(genes_v2) > 0) {
      gene_lists[[paste0("Visit2_", str_to_title(module))]] <- genes_v2
    }
  }
  
  # Visit 3
  if (module %in% unique(visit3_modules)) {
    genes_v3 <- extract_module_genes(visit3_modules, module, Annotation_data)
    if (length(genes_v3) > 0) {
      gene_lists[[paste0("Visit3_", str_to_title(module))]] <- genes_v3
    }
  }
}

# ============================================================================
# STEP 4: Run compareCluster analysis
# ============================================================================

ck_kegg <- compareCluster(
  geneClusters = gene_lists,
  fun = "enrichKEGG",
  organism = "hsa",
  universe = Annotation_data$EntrezGeneID,
  pvalueCutoff = 0.05
)

cat("\ncompareCluster results summary:\n")
print(ck_kegg)

# ============================================================================
# STEP 5: Create FACETED dot plot
# ============================================================================

results_df <- as.data.frame(ck_kegg) %>%
  mutate(
    Visit = str_extract(Cluster, "Visit\\d+"),
    Module = str_extract(Cluster, "(?<=_)[A-Za-z]+$"),
    GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"),
                               function(x) as.numeric(x[1]) / as.numeric(x[2]))
  ) %>%
  filter(!is.na(Visit) & !is.na(Module))

View(results_df)

results_df$Module <- factor(results_df$Module, levels = c("Turquoise", "Blue", "Brown", "Red", "Black"))

# Save results as csv and .rds file
write.csv(results_df, "data/Tracking_pathway_enrichment_results_across_3visits.csv")
saveRDS(results_df, file = "data/results_df_pathway_enrichment.rds")

p_faceted <- ggplot(results_df, aes(x = Visit, y = Description)) +
  geom_point(aes(size = GeneRatio_numeric, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", name = "padj") +
  scale_size_continuous(name = "GeneRatio", range = c(2, 8)) +
  facet_wrap(~ Module, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(x = "Visit", y = "Pathway")

ggsave("data/faceted_dotplot_by_module.png", p_faceted, 
       width = 16, height = 12, dpi = 300)

# ============================================================================
# STEP 6: Top pathways per module (optional)
# ============================================================================

top_pathways_per_module <- results_df %>%
  group_by(Module) %>%
  arrange(p.adjust) %>%
  slice_head(n = 50) %>%
  ungroup()

p_faceted_top <- ggplot(top_pathways_per_module, aes(x = Visit, y = Description)) +
  geom_point(aes(size = GeneRatio_numeric, color = p.adjust)) + 
  scale_color_gradient(low = "red", high = "blue", name = "padj") +
  scale_size_continuous(name = "GeneRatio", range = c(2, 8)) +
  facet_wrap(~ Module, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(x = "Visit", y = "Pathway")


ggsave("data/faceted_dotplot_top50_by_module.png", p_faceted_top, 
       width = 20, height = 24, dpi = 300)


# To extract both geneid and gene name from a module of a visit
extract_pathway_genes <- function(results_df, module, visit, pathway, annotation_data = Annotation_data) {
  entry <- results_df %>%
    filter(Module == module, Visit == visit, Description == pathway)
  
  if (nrow(entry) == 0) return(NULL)
  
  # Split Entrez IDs into vector
  entrez_genes <- strsplit(entry$geneID, "/")[[1]]
  
  # Map to gene symbols using annotation_data
  gene_symbols <- annotation_data %>%
    filter(EntrezGeneID %in% entrez_genes) %>%
    pull(EntrezGeneSymbol)
  
  # Return both Entrez IDs and gene symbols
  return(list(entrez = entrez_genes, symbol = gene_symbols))
}

# Example use - 
Black_Visit1_Complement_and_coagulation_cascades <- extract_pathway_genes(results_df, module = "Black", visit = "Visit1", pathway = "Complement and coagulation cascades")
saveRDS(Black_Visit1_Complement_and_coagulation_cascades, file = "data/Black_Visit1_Complement_and_coagulation_cascades.rds")




#### 3. Quick Cell type enrichment check - EWCE Enrichment (for one module of one visit) - Publication ready ewce enrichment plot code is the next code chunk)----------------------------

bwnet_visit1 <- readRDS("data/visit1_results/bwnet_visit_1.rds")  # load the wgcna network result for a particular visit (as saved in wgcna workflow in WGCNA code script)

module_proteins <- names(bwnet_visit1$colors)[bwnet_visit1$colors == "blue"] # input the module colour for enrichment

Annotation_data <- read_excel("Somalogic_7-5K_protein_list.xlsx") # load annotation code

# Modifying the dataframe to have just SeqId, EntrezGeneID and EntrezGeneSymbol
Annotation_data <- Annotation_data %>%
select("SeqId", "EntrezGeneID", "EntrezGeneSymbol")

# Modifying the seqID column
Annotation_data_new <- Annotation_data %>%
  mutate(SeqId = gsub("^(\\d+)-(\\d+)$", "seq.\\1.\\2", SeqId))

module_proteins <- as.data.frame(module_proteins)
colnames(module_proteins) <- "SeqId"

associated_proteins <- left_join(module_proteins, Annotation_data_new,  by = "SeqId", keep = F)
View(associated_proteins)

ctd <- ewceData::ctd()
hits <- associated_proteins$EntrezGeneSymbol 
reps <- 10000
annotLevel <- 1
## Enrichment
blue_ct <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                           sctSpecies = "mouse",
                                           genelistSpecies = "human",
                                           hits = hits,
                                           reps = reps,
                                           annotLevel = annotLevel)

plot_list <- EWCE::ewce_plot(total_res = blue_ct$results,
                             mtc_method ="BH")
print(plot_list$plain)
### ----




#### 4. Publication ready Cell type enrichment facet plot for modules across visits using EWCE ----------------------------------------------------------------------

# ============================================================================ 
# STEP 1: Load annotation data and module assignments
# ============================================================================

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

# ============================================================================ 
# STEP 2: Function to extract genes from module assignment vector
# ============================================================================

extract_module_genes_vector <- function(module_vector, module_color, annotation_data) {
  module_genes <- names(module_vector)[module_vector == module_color]
  module_genes_df <- data.frame(SeqId = module_genes, stringsAsFactors = FALSE)
  associated_genes <- left_join(module_genes_df, annotation_data, by = "SeqId")
  entrez_genes <- na.omit(associated_genes$EntrezGeneSymbol)
  return(entrez_genes)
}

# ============================================================================ 
# STEP 3: Function to run EWCE per module
# ============================================================================

perform_ewce_vector <- function(module_vector, module_color, annotation_data, ctd_data, visit_name, reps = 10000, annotLevel = 1) {
  genes <- extract_module_genes_vector(module_vector, module_color, annotation_data)
  if(length(genes) < 5) return(NULL)
  
  res <- tryCatch(
    EWCE::bootstrap_enrichment_test(
      sct_data = ctd_data,
      sctSpecies = "mouse",
      genelistSpecies = "human",
      hits = genes,
      reps = reps,
      annotLevel = annotLevel
    ),
    error = function(e) NULL
  )
  
  if(is.null(res)) return(NULL)
  
  df <- res$results
  df$Module <- module_color
  df$Visit <- visit_name
  df$n_genes <- length(genes)
  return(df)
}

# ============================================================================ 
# STEP 4: Run EWCE for all visits/modules using module assignment vectors
# ============================================================================

run_ewce_all_visits <- function(module_vectors_list, annotation_data, ctd_data, reps = 10000, annotLevel = 1) {
  all_results <- list()
  
  for(visit_name in names(module_vectors_list)) {
    message("Processing ", visit_name)
    mod_vec <- module_vectors_list[[visit_name]]
    for(mod_color in unique(mod_vec)) {
      if(mod_color == "grey") next
      res <- perform_ewce_vector(mod_vec, mod_color, annotation_data, ctd_data, visit_name, reps, annotLevel)
      if(!is.null(res)) all_results[[paste(visit_name, mod_color, sep = "_")]] <- res
    }
  }
  
  combined <- bind_rows(all_results)
  return(combined)
}

# ============================================================================ 
# STEP 5: Prepare module assignment vectors and run
# ============================================================================

visit1_modules <- bwnet_visit1_modassignments
visit2_modules <- visit2_modassignments_matched
visit3_modules <- visit3_modassignments_matched_ref_visit2_modassignments_matched


visit_data_list <- list(
  "Visit1" = visit1_modules,
  "Visit2" = visit2_modules,
  "Visit3" = visit3_modules
)

ctd <- ewceData::ctd() # load reference single-cell transcriptome

ewce_combined_results <- run_ewce_all_visits(visit_data_list, Annotation_data, ctd)
View(ewce_combined_results)
saveRDS(ewce_combined_results, file = "data/ewce_combined_results.rds")


# Check if EWCE's q matches manual BH correction
ewce_combined_results <- ewce_combined_results %>%
  mutate(
    q_manual_within = p.adjust(p, method = "BH"), # BH across all rows
    q_manual_by_module = ave(p, interaction(Visit, Module), FUN = function(x) p.adjust(x, method = "BH"))
  )

# Compare the built-in q with your recalculated versions
head(
  ewce_combined_results %>% 
    select(Visit, Module, CellType, p, q, q_manual_within, q_manual_by_module)
)

# ============================================================================ 
# STEP 6: Convert to facet plot format (like pathway enrichment)
# ============================================================================

plot_df <- ewce_combined_results %>%
  mutate(
    Visit = factor(Visit, levels = c("Visit1","Visit2","Visit3")),
    Module = factor(str_to_title(Module),
                     levels = c("Turquoise","Blue","Brown","Black","Yellow")),
    CellType = str_replace_all(CellType, "_", " ")
  )%>%
  filter(q < 0.05)  # only significant cell types

# Example: facet per Module, x = Visit, y = CellType
ewce_ggplot <- ggplot(plot_df, aes(x = Visit, y = CellType)) +
  geom_point(aes(size = -log10(q), color = fold_change)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(2,8)) +
  facet_wrap(~ Module, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 8) +
  theme(

  axis.title.x = element_text(face = "bold", size = 10),
  axis.title.y = element_text(face = "bold", size = 10),
  axis.text.y = element_text(size = 8),
  strip.text = element_text(face = "bold", size = 10),
  legend.title = element_text(face = "bold", size = 9),
  # axis.text.y = element_text(size = 8),
  panel.spacing.y = unit(0.02, "lines")
  ) +
  labs(x = "Visit", y = "Cell Type", color = "Fold Change", size = "-log10(q)")

ggsave("data/ewce_plot.png", ewce_ggplot, width = 8, height = 4, dpi = 300)




#### 5. Extracting proteins of each module of a visit (one visit at a time - repeat entire code chunk for other visits) ---------------------------------------------------
# Load expression data 
counts_data_visit3 <- readRDS("data/visit3_results/counts_data_visit_3.rds") # samples x proteins matrix/data.frame used in WGCNA
# Load module assignments
visit3_modassignments_matched_ref_visit2_modassignments_matched <- 
readRDS("data/visit3_results/visit3_modassignments_matched_ref_visit2_modassignments_matched.rds")

moduleColors <- visit3_modassignments_matched_ref_visit2_modassignments_matched # module assignments of a visit
stopifnot(length(moduleColors) == ncol(counts_data_visit3))  # sanity check

# Load annotation data 
Annotation_data <- read_excel("Somalogic_7-5K_protein_list.xlsx") %>%
  select("SeqId", "EntrezGeneSymbol","EntrezGeneID") %>%
  mutate(SeqId = gsub("^(\\d+)-(\\d+)$", "seq.\\1.\\2", SeqId))

# table of sizes (optional)
sort(table(moduleColors), decreasing = TRUE)

# per-protein table
protein_by_module <- data.frame(
  protein = colnames(counts_data_visit3),
  module  = moduleColors,
  stringsAsFactors = FALSE
)

# # split into a named list; drop "grey" (unassigned) if you want
proteins_in_module <- split(protein_by_module$protein, protein_by_module$module)
proteins_in_module <- proteins_in_module[setdiff(names(proteins_in_module), "grey")]

# Save to .rds 
saveRDS(proteins_in_module, file = "data/proteins_in_relabelled_modules_visit_3.rds")

# example: proteins in the blue module
blue_proteins <- proteins_in_module[["blue"]]

# Have gene names instead of just proteins
proteins_in_modules_visit_3 <- readRDS("data/proteins_in_relabelled_modules_visit_3.rds")

names(proteins_in_modules_visit_3)
lengths(proteins_in_modules_visit_3)

# Make lookup vectors
lookup_symbol <- setNames(Annotation_data$EntrezGeneSymbol, Annotation_data$SeqId)
lookup_entrez <- setNames(Annotation_data$EntrezGeneID, Annotation_data$SeqId)

# Map both GeneSymbol and Entrez IDs for each module
proteins_in_modules_visit_3_annotated <- lapply(proteins_in_modules_visit_3, function(seq_ids) {
  df <- data.frame(
    SeqId = seq_ids,
    GeneSymbol = lookup_symbol[seq_ids],
    EntrezGeneID = lookup_entrez[seq_ids],
    stringsAsFactors = FALSE
  )
  rownames(df) <- NULL
  df
})

saveRDS(proteins_in_modules_visit_3_annotated, file = "data/proteins_in_relabelled_modules_visit_3_annotated.rds")

### ----


#### 6. TISSUE ENRICHMENT using Enrichr for one module of one visit ----------------------------------------------------------------------
# Load the gene names of a module of any visit (this can be generated using the code in code chunk 5)
genes_visit1_modules <- readRDS("data/proteins_in_modules_visit_1_annotated.rds")
View(genes_visit1_modules)
my_genes <- genes_visit1_modules$green$GeneSymbol # choose the module 

dbs <- c("Human_Gene_Atlas", "ARCHS4_Tissues","Allen_Brain_Atlas_10x_scRNA_2021", 
"Allen_Brain_Atlas_Up","Allen_Brain_Atlas_Down")
res_er <- enrichr(my_genes, dbs)  # plain data.frame output

# Summary function
summarize_tissue_enrichment <- function(res) {
  for(db in names(res)) {
    cat("\n==========", db, "==========\n")
    sig <- res[[db]][res[[db]]$Adjusted.P.value < 0.05, ]
    if(nrow(sig) > 0) {
      top20 <- head(sig[, c("Term", "Overlap", "P.value","Adjusted.P.value", "Combined.Score")], 100) #, "Genes"
      print(top20)
    } else {
      cat("No significant enrichments\n")
    }
  }
}

summarize_tissue_enrichment(res_er)
### ----

#### 7. TISSUE ENRICHMENT using Enrichr for all the modules of all the visits ------

# ============================================================================
# STEP 1: Load pre-extracted gene modules for all visits
# ============================================================================

genes_visit1_modules <- readRDS("data/proteins_in_modules_visit_1_annotated.rds")
genes_visit2_modules <- readRDS("data/proteins_in_relabelled_modules_visit_2_annotated.rds")
genes_visit3_modules <- readRDS("data/proteins_in_relabelled_modules_visit_3_annotated.rds")

# ============================================================================
# STEP 2: Define tissue databases
# ============================================================================

tissue_dbs <- c( "ARCHS4_Tissues","Allen_Brain_Atlas_10x_scRNA_2021", 
"Allen_Brain_Atlas_Up","Allen_Brain_Atlas_Down")

# ============================================================================
# STEP 3: Define modules of interest
# ============================================================================

modules_of_interest <- c("turquoise", "blue", "brown", "red", "black", "yellow", "green")

# ============================================================================
# STEP 4: Run EnrichR for all modules and visits
# ============================================================================

all_results <- list()

# Visit 1
cat("Processing Visit 1...\n")
for (module in modules_of_interest) {
  if (module %in% names(genes_visit1_modules)) {
    my_genes <- genes_visit1_modules[[module]]$GeneSymbol
    
    if (length(my_genes) > 10) {  # Need sufficient genes for enrichment
      cat("  Running EnrichR for", module, "module (", length(my_genes), "genes)\n")
      res_er <- enrichr(my_genes, tissue_dbs)
      
      # Combine all database results with metadata
      for (db_name in names(res_er)) {
        df <- res_er[[db_name]]
        if (nrow(df) > 0) {
          df$Visit <- "Visit1"
          df$Module <- str_to_title(module)
          df$Database <- db_name
          df$Cluster <- paste0("Visit1_", str_to_title(module))
          all_results[[paste0("Visit1_", module, "_", db_name)]] <- df
        }
      }
    }
  }
}

# Visit 2
cat("\nProcessing Visit 2...\n")
for (module in modules_of_interest) {
  if (module %in% names(genes_visit2_modules)) {
    my_genes <- genes_visit2_modules[[module]]$GeneSymbol
    
    if (length(my_genes) > 10) {
      cat("  Running EnrichR for", module, "module (", length(my_genes), "genes)\n")
      res_er <- enrichr(my_genes, tissue_dbs)
      
      for (db_name in names(res_er)) {
        df <- res_er[[db_name]]
        if (nrow(df) > 0) {
          df$Visit <- "Visit2"
          df$Module <- str_to_title(module)
          df$Database <- db_name
          df$Cluster <- paste0("Visit2_", str_to_title(module))
          all_results[[paste0("Visit2_", module, "_", db_name)]] <- df
        }
      }
    }
  }
}

# Visit 3
cat("\nProcessing Visit 3...\n")
for (module in modules_of_interest) {
  if (module %in% names(genes_visit3_modules)) {
    my_genes <- genes_visit3_modules[[module]]$GeneSymbol
    
    if (length(my_genes) > 10) {
      cat("  Running EnrichR for", module, "module (", length(my_genes), "genes)\n")
      res_er <- enrichr(my_genes, tissue_dbs)
      
      for (db_name in names(res_er)) {
        df <- res_er[[db_name]]
        if (nrow(df) > 0) {
          df$Visit <- "Visit3"
          df$Module <- str_to_title(module)
          df$Database <- db_name
          df$Cluster <- paste0("Visit3_", str_to_title(module))
          all_results[[paste0("Visit3_", module, "_", db_name)]] <- df
        }
      }
    }
  }
}

# ============================================================================
# STEP 5: Combine all results and prepare for plotting
# ============================================================================

combined_results <- bind_rows(all_results) %>%
  filter(Adjusted.P.value < 0.05) %>%
  mutate(
    # Convert Overlap to GeneRatio (proportion) for better visualization
    # Overlap format: "5/200" -> GeneRatio: 0.025
    Overlap_percentage = sapply(strsplit(as.character(Overlap), "/"),
                               function(x) as.numeric(x[1]) / as.numeric(x[2])),
    # Extract just the count for alternative visualization
    Overlap_count = sapply(strsplit(as.character(Overlap), "/"),
                          function(x) as.numeric(x[1])),
    Module = factor(Module, levels = c("Turquoise", "Blue", "Brown", "Red", "Black", "Yellow", "Green"))
  ) %>%
  rename(
    Description = Term,
    p.adjust = Adjusted.P.value
  )

# Save combined results
write.csv(combined_results, "data/Tissue_enrichment_results_across_3visits.csv", row.names = FALSE)
saveRDS(combined_results, "data/results_df_tissue_enrichment.rds")

cat("\n=== PROCESSING COMPLETE ===\n")
cat("Total significant tissue enrichments:", nrow(combined_results), "\n\n")

# View sample of results
cat("Preview of results:\n")
print(head(combined_results %>% select(Visit, Module, Database, Description, p.adjust, Overlap_percentage)))

# ============================================================================
# STEP 6: Create faceted dotplot (all significant results)
# ============================================================================

# Filter out mouse related enrichments
combined_results <- combined_results %>% 
filter(!grepl("(?i)mouse|mus musculus", Description))


if (nrow(combined_results) > 0) {
  p_faceted <- ggplot(combined_results, aes(x = Visit, y = Description)) +
    geom_point(aes(size = Overlap_percentage, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "padj") +
    scale_size_continuous(name = "Overlap_percentage", range = c(2, 8)) +
    facet_wrap(~ Module, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    labs(x = "Visit", y = "Cell Type or Tissue") #title = "Tissue/Organ Enrichment Across Modules and Visits")
  
  ggsave("data/faceted_dotplot_tissue_enrichment_by_module.png", 
         p_faceted, width = 20, height = 24, dpi = 300)
  
  cat("\nPlot saved: faceted_dotplot_tissue_enrichment_by_module.png\n")
}

# ============================================================================
# STEP 7: Top tissues per module (recommended for cleaner visualization)
# ============================================================================

if (nrow(combined_results) > 0) {
  top_tissues_per_module <- combined_results %>%
    group_by(Module) %>%
    arrange(p.adjust) %>%
    slice_head(n = 30) %>%  # Top 30 per module
    ungroup()
  
  p_faceted_top <- ggplot(top_tissues_per_module, aes(x = Visit, y = Description)) +
    geom_point(aes(size = Overlap_percentage, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "padj") +
    scale_size_continuous(name = "GeneRatio", range = c(2, 8)) +
    facet_wrap(~ Module, scales = "free_y", ncol = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    labs(x = "Visit", y = "Tissue/Organ", title = "Top 30 Tissue Enrichments per Module")
  
  ggsave("data/faceted_dotplot_top30_tissue_enrichment_by_module.png", 
         p_faceted_top, width = 20, height = 20, dpi = 300)
  
  cat("Plot saved: faceted_dotplot_top30_tissue_enrichment_by_module.png\n")
}



#### 8. Protein transitions between modules across visits (Sankey/Alluvial plots) ---------------------------------------------

# Load module assignments of each visit
bwnet_visit1 <- readRDS("data/visit1_results/bwnet_visit_1.rds")
bwnet_visit1_modassignments <- bwnet_visit1$colors

visit2_modassignments_matched <- readRDS("data/visit2_results/visit2_modassignments_matched.rds")

visit3_modassignments_matched_ref_visit2_modassignments_matched <- readRDS("data/visit3_results/visit3_modassignments_matched_ref_visit2_modassignments_matched.rds")

# Create a data frame called module_assignments where each row represents a protein and its corresponding module assignments across three visits.
module_assignments <- data.frame(
  Protein = names(bwnet_visit1_modassignments),
  Visit1 = bwnet_visit1_modassignments,
  Visit2 = visit2_modassignments_matched,
  Visit3 = visit3_modassignments_matched_ref_visit2_modassignments_matched
)

head(module_assignments)

# Create a summary of transitions between module assignments across the three visits.
transitions <- module_assignments %>%
  group_by(Visit1, Visit2, Visit3) %>% # Group the data by combinations of module assignments across the three visits.
  summarise(Count = n(), .groups = 'drop') # For each combination of Visit1, Visit2, and Visit3, count the number of occurrences and store it in a column called Count. The .groups = 'drop' argument ensures that the grouping structure is dropped after summarization.

head(transitions)


# Build stratum data separately
p_base <- ggplot(transitions,
                 aes(axis1 = Visit1, axis2 = Visit2, axis3 = Visit3, y = Count)) +
  geom_alluvium(aes(fill = Visit1), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "grey")

strata_data <- ggplot_build(p_base)$data[[2]]  

# Split by axis
visit1_labels <- subset(strata_data, x == 1)   # axis1
visit2_labels <- subset(strata_data, x == 2)   # axis2
visit3_labels <- subset(strata_data, x == 3)   # axis3

# Final plot
p_final <- p_base +
  # Left-side labels for Visit1
  geom_label_repel(
    data = visit1_labels,
    aes(x = x, y = y, label = stratum),
    nudge_x = -0.2, direction = "y",
    size = 4.2, segment.color = "grey40",
    inherit.aes = FALSE
  ) +
  # Right-side labels for Visit3
  geom_label_repel(
    data = visit3_labels,
    aes(x = x, y = y, label = stratum),
    nudge_x = 0.2, direction = "y",
    size = 4.2, segment.color = "grey40",
    inherit.aes = FALSE
  ) +
  # Smaller inline labels for Visit2
  geom_text(
    data = visit2_labels,
    aes(x = x, y = y, label = stratum),
    size = 3.2, color = "black",
    inherit.aes = FALSE
  ) +
  scale_x_discrete(limits = c("Visit1", "Visit2", "Visit3"),
                   expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  labs(
       y = "Number of Proteins") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# Save
png("data/Protein_Module_Transitions_Across_Visits_Sankey.png",
    width = 3000, height = 4000, res = 300)
print(p_final)
dev.off()
