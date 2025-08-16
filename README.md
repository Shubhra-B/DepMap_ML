# Dependency_Map_wrappers
Depmap has omics dataset from cell lines modelling various malignancies. 

There are a few online/point-and-click applications but have limited customisation options to visualise according to needs and EDA (see https://pmc.ncbi.nlm.nih.gov/articles/PMC7924953/). 

Here, there are functions which can be readily used  irectly on rawdata allowing customisation in terms of not only visualisation but also statistical analysis**.

**TPM_visualisation.R** : TPM normalised transcripts for specific genes to visually inspect difference between  cell lines.

**Filter_visualise_mutations.R** : Filter mutations in specific cell lines based on predicted damaging effect and plot on lollipop plot.

**calculate_rotations_pca.R** : Perform PCA for specific group of cell lines with specific mutations (you may use the scripts above/or any other criteria) and shortlist 50 top genes                                     which are contributing to the variability of PC1 and PC2.
                                Plots PCA plot as well.


