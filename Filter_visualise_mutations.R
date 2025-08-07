
library(GenomicRanges)
library(tidyverse)
library(dandelion)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(trackViewer)
library(RColorBrewer)


#filter only damaging mutations, it is not uncommon that a single cell line has more than one mutation.

#use the example file to get the appropriate column names

filter_deleterious_mutations <- function(df, gene) {
  df_filtered <- df[apply(df, 1, function(row) {
    df<-df %>% select(ModelID,HugoSymbol, HugoSymbol, LikelyLoF,Sift,Polyphen,VepImpact,VariantInfo,Chrom,Pos,Ref,Alt,AF,DP)
    any(grepl("deleterious|damaging|high", row, ignore.case = TRUE)) &&
      any(grepl(gene, row, ignore.case = TRUE))
  }), ]
  
  return(df_filtered)
}

#example usage
STAG2_damaging <- filter_deleterious_plus_gene(somatic_mut_cleaned, "STAG2")


#Plot the mutations, we use dandelion plot since it is more appropriate for densely distributed mutations

plot_mutations <- function(gene_name, df,chromosome, output_dir = ".") {
  gr2 <- GRanges(
    chromosome, 
    IRanges(
      min(df$Pos) - 500, 
      max(df$Pos) + 500, 
      names = gene_name))

  #SNP data needs to be gRanges
  SNP_data <- df %>% dplyr::select(Chrom, Pos, ModelID)
  SNP_data$Pos2 <- SNP_data$Pos
  SNP_data <- SNP_data[, c(1, 2, 4, 3)]
  colnames(SNP_data) <- c("chrom", "chromStart", "chromEnd", "name")
  SNP_data_granges <- GenomicRanges::makeGRangesFromDataFrame(SNP_data)

  #transcript model: Make sure you know if you want hg19 or grch38
  transcript_model <- geneModelFromTxdb(
    TxDb.Hsapiens.UCSC.hg38.knownGene,
    org.Hs.eg.db,
    gr = gr2
  )

  # Extract features
  features_to_plot <- c(range(transcript_model[[1]]$dat), range(transcript_model[[2]]$dat))
  names(features_to_plot) <- c(transcript_model[[1]]$name, transcript_model[[2]]$name)

  file_name <- paste0(gene_name, "_dandelion.pdf")
  plot_path <- file.path(output_dir, file_name)

  pdf(plot_path, width = 10, height = 8)
  dandelion.plot(SNP_data_granges, features_to_plot, ranges = gr2, type = "pin")
  dev.off()

}

#example usage
plot_mutations("STAG2", STAG2_damaging, "chrX")

##########If you would like to see 
#1. Which are the lineage of cancer cell lines carrying these damaging mutation in your gene of interest
#2. What percentage of cell lines from your lineage of interest  
#3. Within your gene of interest and lineage of interest what are the disease which are being modelled.


plot_lineage_damaging_percentage <- function(mutation_df, model_path = "../Model.csv", lineage_of_interest) {
mut_IDs <- unique(mutation_df$ModelID)
model_cleaned <- read.csv(model_path, header = TRUE) %>%
select(ModelID, OncotreeLineage, OncotreePrimaryDisease)

summary_merged <- model_cleaned %>%
mutate(damaging = ModelID %in% mut_IDs) %>%
group_by(OncotreeLineage) %>%
                          summarise(
                          Cell_line_with_damaging = sum(damaging),
                          Total_number_of_celllines_in_Lineage = n(),
                          .groups = "drop") %>%

mutate(percentage_with_damaging = 100 * Cell_line_with_damaging / Total_number_of_celllines_in_Lineage)
  
p <- ggplot(summary_merged, aes(x = reorder(OncotreeLineage, -percentage_with_damaging), y = percentage_with_damaging)) +
geom_col() +
theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_text(aes(label = Total_number_of_celllines_in_Lineage), vjust = -0.25)
  
    print(p)
    lineage_disease_counts <- model_cleaned %>%
    filter(OncotreeLineage == lineage_of_interest) %>%
    count(OncotreePrimaryDisease)
  
    # Plot disease counts for the lineage_of_interest
    p_disease <- ggplot(lineage_disease_counts, aes(x = OncotreePrimaryDisease, y = n)) +
    theme_classic() +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Distribution of OncotreePrimaryDisease in", lineage_of_interest),
         x = "Oncotree Primary Disease",
         y = "Count")
  
  print(p_disease)
  
  return(list(
    summary = summary_merged,
    damaging_plot = p,
    lineage_disease_counts = lineage_disease_counts,
    disease_distribution_plot = p_disease
  ))
}
#this give us all the cellines with STAG2 mutations in myeloid lineage and disease which are present in the cell lines of lineage of interest
result <- plot_lineage_damaging_percentage(STAG2_damaging, lineage_of_interest = "Myeloid")
result$damaging_plot         # damaging mutation percentage plot
result$disease_distribution_plot  # OncotreePrimaryDisease plot


