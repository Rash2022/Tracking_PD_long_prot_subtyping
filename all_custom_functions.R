# ALL CUSTOM FUNCTIONS


# 1. FOR GRABBING EXPRESSION DATA/APTAMERS ONLY 
# Used in wgcna code script
is_seq <- function(.x) grepl("^seq\\.[0-9]{4}", .x) # regex for APTAMERS


# 2. FOR AUTOMATING PATHWAY ENRICHMENT FOR ALL MODULES OF ANY VISIT IN PRE-DEFINED OUTPUT FOLDER 
# Used in the enrichment code script

output_dir <- "data/results"

perform_pathway_analysis <- function(module_color, bwnet, annotation_data, output_prefix) {
  # Extract proteins for the given module color 
  module_indices <- which(bwnet$colors == module_color)
  
  # Check if any proteins in this module
  if(length(module_indices) == 0) {
    message(paste("No proteins found for module color:", module_color))
    return(NULL)
  }
  
  # Get protein names from the indices
  module_proteins <- data.frame(SeqId = names(bwnet$colors)[module_indices], 
                                stringsAsFactors = FALSE)
  
  # Join with annotation data to get EntrezGeneIDs
  associated_proteins <- left_join(module_proteins, annotation_data, by = "SeqId", keep = FALSE)
  
  # Skip if no or too few proteins with EntrezGeneIDs
  valid_genes <- na.omit(associated_proteins$EntrezGeneID)
  if(length(valid_genes) < 5) {
    message(paste("Module", module_color, "has fewer than 5 proteins with EntrezGeneIDs. Skipping..."))
    return(NULL)
  }
  
  # Perform KEGG enrichment
  tryCatch({
    kegg_result <- enrichKEGG(
      gene = valid_genes,
      organism = "hsa",
      universe = annotation_data$EntrezGeneID,
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "fdr",
      use_internal_data = TRUE
    )
    
    # Check if any enriched pathways were found
    if(nrow(kegg_result@result) > 0) {
      # Define output folder
      output_dir <- output_dir
      # Save KEGG barplot inside the folder
      # output_file <- paste0(output_prefix, "_KEGG_barplot_", module_color, ".png")
      output_file <- file.path(output_dir, paste0(output_prefix, "_KEGG_barplot_", module_color, ".png"))
      png(output_file, width = 1200, height = 2000, res = 150)
      print(barplot(kegg_result, showCategory = min(40, nrow(kegg_result@result))))
      dev.off()
      message(paste("KEGG analysis for module", module_color, "saved to", output_file))
      
      # Return the result for potential further analysis
      return(kegg_result)
    } else {
      message(paste("No significant KEGG pathways found for module", module_color))
      return(NULL)
    }
  }, error = function(e) {
    message(paste("Error in KEGG analysis for module", module_color, ":", e$message))
    return(NULL)
  })
}

# Main function to analyze all modules
analyze_all_modules <- function(bwnet, annotation_file, output_prefix) {
  # Load and process annotation data
  Annotation_data <- read_excel(annotation_file)
  Annotation_data_new <- Annotation_data %>%
    select("SeqId", "EntrezGeneID") %>%
    mutate(SeqId = gsub("^(\\d+)-(\\d+)$", "seq.\\1.\\2", SeqId))
  
  # Get all module colors
  module_colors <- unique(bwnet$colors)
  
  # Add debugging information
  message("Module colors found: ", paste(module_colors, collapse=", "))
  message("Structure of bwnet$colors: ")
  print(str(bwnet$colors))
  
  # Create a list to store results
  all_results <- list()
  
  # Process each module
  for(color in module_colors) {
    message(paste("Processing module:", color))
    result <- perform_pathway_analysis(color, bwnet, Annotation_data_new, output_prefix)
    if(!is.null(result)) {
      all_results[[color]] <- result
    }
  }
  
  # Return all results
  return(all_results)
}
