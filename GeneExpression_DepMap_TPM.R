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


############Funciton that takes gene names and filters only tpm of the user defined genes####################

genes_of_interest=c("X","STAG1","STAG2","RAD21","SMC1A","SMC3","WAPL","PDS5A","PDS5B") #change here

tpm_selection <- function(df, genes) {
tpm1=tpm[,colnames(df) %in% genes, ]
colname_tpm1 <- paste0(colnames(tpm1), "_RNA")
colnames(tpm1)<-colname_tpm1
return(tpm1)
 }

#example usage
TPM_cohesin_components <- tpm_selection(tpm, genes_of_interest)
TPM_cohesin_components<-column_to_rownames(TPM_cohesin_components,"X_RNA")


tpm_violin<-ggplot(df_TPM_cohesin_components, aes(gene, TPM)) + geom_violin()

