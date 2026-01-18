#' Run Tissue Enrichment Analysis
#'
#' @param input_genes Vector of gene symbols.
#' @param background_genes Vector of background gene symbols.
#' @param specificity_types Character vector of HPA categories to use.
#' @export
run_tissue_enrichment <- function(input_genes, 
                                  background_genes = NULL, 
                                  specificity_types = c("Tissue enriched", "Group enriched", "Tissue enhanced")) {
  
  # 1. Load Internal Data
  if (!exists("hpa_dataset")) {
    stop("Package data 'hpa_dataset' not found.")
  }
  
  # 2. Setup Inputs
  input_genes <- unique(toupper(input_genes))
  
  if (is.null(background_genes)) {
    # If no background, fallback to all HPA genes (same as your script)
    background_genes <- unique(hpa_dataset$Gene)
  } else {
    background_genes <- unique(toupper(background_genes))
  }

  # 3. Filter Reference Data
  # Only keep HPA entries that match the requested categories
  universe_map <- hpa_dataset %>%
    filter(Specificity %in% specificity_types) %>%
    # Only keep HPA entries that exist in your background
    filter(Gene %in% background_genes)
  
  # --- LOGIC SYNC WITH YOUR SCRIPT ---
  
  # N: Total size of the background universe
  # Your script uses the size of the *provided background*, not the intersection.
  N <- length(background_genes)
  
  # K: Total size of the target list found in background
  # Your script intersects input with background (not HPA)
  valid_input_genes <- intersect(input_genes, background_genes)
  K <- length(valid_input_genes)
  
  # -----------------------------------
  
  if (K == 0) stop("No overlap between input genes and background genes.")

  # 4. Run Stats per Tissue
  all_tissues <- unique(universe_map$Tissue)
  results_list <- list()
  
  for (tissue in all_tissues) {
    
    # M: Number of specific genes for this tissue in the HPA (filtered by background)
    tissue_specific_genes <- universe_map %>%
      filter(Tissue == tissue) %>%
      pull(Gene) %>%
      unique()
    
    M <- length(tissue_specific_genes)
    
    # x: Overlap (Input genes that are Tissue Specific)
    overlap_genes <- intersect(valid_input_genes, tissue_specific_genes)
    x <- length(overlap_genes)
    
    if (x > 0) {
      # Fold Change: (x/K) / (M/N)
      fold_change <- (x / K) / (M / N)
      
      # P-value: Hypergeometric
      p_val <- phyper(x - 1, M, N - M, K, lower.tail = FALSE)
      
      results_list[[tissue]] <- data.frame(
        Tissue = tissue,
        Tissue_Specific_Genes_Count = M,
        Overlap_Count = x,
        Fold_Change = fold_change,
        PValue = p_val,
        Genes = paste(overlap_genes, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # 5. Final Output
  if (length(results_list) == 0) return(NULL)
  
  final_df <- do.call(rbind, results_list)
  final_df$Adj_PValue <- p.adjust(final_df$PValue, method = "BH")
  
  # --- IMPROVEMENT: Reorder columns so Genes is last and Adj_P is near P ---
  final_df <- final_df %>%
    select(
      Tissue, 
      Tissue_Specific_Genes_Count, 
      Overlap_Count, 
      Fold_Change, 
      PValue, 
      Adj_PValue,   # Moved here
      Genes         # Moved to the end
    ) %>%
    arrange(Adj_PValue) # Sort by significance
  
  return(final_df)
}