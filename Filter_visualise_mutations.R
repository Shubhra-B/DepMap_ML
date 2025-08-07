
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


