# Load required libraries
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(readr)
library(GO.db)

# ---- INPUT FILES ----
# Get all CSV files in report folder that contain "vs" in the name
report_files <- list.files("reports", pattern = "_vs_.*\\.csv$", full.names = TRUE)
report_files
# Read GO annotation from GAF file (columns: gene_id in col 2, GO ID in col 5)
gaf_file <- "lib/FungiDB-65_AfumigatusA1163_GO.gaf"
go_table <- read_tsv(gaf_file, comment = "!", col_names = FALSE) %>%
    select(term = X5, gene = X2)

# Get all GO IDs and terms as a named list/array
go_data <- as.list(GOTERM)

# Extract GO IDs (names of the list) and terms (values of the list)
go_ids <- names(go_data)
go_terms <- sapply(go_data, Definition)

# Create a data frame
go_df <- data.frame(
  GO_ID = go_ids,
  GO_Term = go_terms,
  stringsAsFactors = FALSE
)

# Write to a tab-delimited file
#write.table(go_df, file = "GO_ID_to_Term.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#t="reports/Azole_genotype_OE_vs_Comp.csv"
#genes_table <- read_csv(t) %>% mutate(gene_id = gsub('-T$','',`...1` ))

# ---- PROCESS EACH FILE ----
for (genes_file in report_files) {
    # Extract output prefix by removing .csv extension
    output_prefix <- sub("\\.csv$", "", basename(genes_file))
    
    cat("Processing:", genes_file, "\n")
    
    # genes_table: columns = gene_id, logFC, pvalue, padj
    genes_table <- read_csv(genes_file) %>% mutate(gene_id = gsub('-T$','',`...1` ))
    
    # ---- PREPARE GENE LISTS ----
    # Define up/down regulated genes (adjust thresholds as needed)
    up_genes <- genes_table %>% filter(log2FoldChange > 1, padj < 0.05) %>% pull(gene_id)
    down_genes <- genes_table %>% filter(log2FoldChange < -1, padj < 0.05) %>% pull(gene_id)
    all_genes <- genes_table$gene_id
    
    # Skip if no significant genes
    if (length(up_genes) == 0 && length(down_genes) == 0) {
        cat("  No significant genes found, skipping...\n")
        next
    }
    
    # ---- ENRICHMENT ANALYSIS ----
    # Use enricher from clusterProfiler with custom gene2GO mapping
    enrich_up <- enricher(gene=up_genes, TERM2GENE = go_table,  
                          #TERM2NAME = go_df,
                          universe = all_genes, 
                          pvalueCutoff  = 0.05,
                          pAdjustMethod = "BH")
    enrich_down <- enricher(gene=down_genes, TERM2GENE = go_table, 
                            #TERM2NAME = go_df,
                            universe = all_genes,
                            pAdjustMethod = "BH"
                            )
    
    # ---- PLOT RESULTS ----
    # Combine and label
    up_df <- as.data.frame(enrich_up@result) %>% mutate(Regulation = "Up")
    down_df <- as.data.frame(enrich_down@result) %>% mutate(Regulation = "Down")
    plot_df <- bind_rows(up_df, down_df) %>%
        filter(p.adjust < 0.05) %>%
        arrange(p.adjust) %>%
        group_by(Regulation) %>%
        slice_head(n = 10) # Top 10 per group
    
    # Skip if no enriched terms
    if (nrow(plot_df) == 0) {
        cat("  No enriched GO terms found, skipping...\n")
        next
    }
    
    # Plot
    p <- ggplot(plot_df, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust), fill = Regulation)) +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() +
        labs(x = "GO Term", y = "-log10(adj. p-value)", title = "Differentially Enriched GO Terms") +
        theme_bw()
    
    # ---- SAVE PLOT ----
    ggsave(paste0(output_prefix, "_differential_GO_terms.png"), plot = p, width = 24, height = 6)
    ggsave(paste0(output_prefix, "_differential_GO_terms.pdf"), plot = p, width = 24, height = 6)
    
    cat("  Saved:", paste0(output_prefix, "_differential_GO_terms.png/pdf\n"))
}
