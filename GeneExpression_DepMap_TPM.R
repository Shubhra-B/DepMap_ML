library("tidymodels")
library("readr")
library("broom.mixed")
library("dotwhisker")
library("recipes")
library("modeldata")
library("themis")
library("RColorBrewer")

############Get data####################

wget "https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsExpressionProteinCodingGenesTPMLogp1.csv"
wget "https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q2&filename=OmicsSomaticMutations.csv"

tpm=read.csv("OmicsExpressionProteinCodingGenesTPMLogp1.csv",header=T)
somatic_mut=read.csv("OmicsSomaticMutations.csv",header=T) 
colnames(tpm)=sub("\\..*","",colnames(tpm))

tpm=read.csv("/home/sbhattacharya/Documents/7.DepMap/OmicsExpressionProteinCodingGenesTPMLogp1.csv",header=T)
somatic_mut=read.csv("/home/sbhattacharya/Documents/7.DepMap/OmicsSomaticMutations.csv",header=T) 
colnames(tpm)=sub("\\..*","",colnames(tpm))
tpm<-column_to_rownames(tpm,"X")

############Function that takes gene names and filters only tpm of the user defined genes####################

genes_of_interest=c("STAG1","STAG2","RAD21","SMC1A","SMC3","WAPL","PDS5A","PDS5B") #change here

tpm_selection <- function(df, genes) {
tpm1=tpm[,colnames(df) %in% genes, ]
return(tpm1)
 }

#example usage
TPM_cohesin_components <- tpm_selection(tpm, genes_of_interest)

func_violin_dist <- function(df) {
  df <- rownames_to_column(df, "CellID")
  df <- df %>%
  pivot_longer(!CellID, names_to = "Genes", values_to = "TPM")
  violin <- ggplot(df, aes(Genes, TPM, fill = Genes)) + 
  geom_violin()
  return(violin)
}
func_violin_dist(TPM_cohesin_components)
