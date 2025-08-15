get_top_pca_genes <- function(data, pc = "PC1", n = 50, file = "pca_top50.txt") {
  pca <- prcomp(data, center = TRUE)
  top <- sort(abs(pca$rotation[, pc]), TRUE)[1:n]
  out <- data.frame(gene = names(top), rotation = top)
  write.table(out, file, sep = "\t", row.names = FALSE, quote = FALSE)
  out}
get_top_pca_genes(tpm_only_myeloid)
run_myeloid_pca <- function(data, coldata, pc = "PC1", n = 50, file = "pca_top50.txt") {
  pca <- prcomp(data, center = TRUE)
  top <- sort(abs(pca$rotation[, pc]), TRUE)[1:n]
  write.table(data.frame(gene = names(top), rotation = top), file, sep = "\t", row.names = FALSE, quote = FALSE)
  eig <- tibble(PC = factor(seq_along(pca$sdev)), variance = pca$sdev^2) %>% mutate(pct = variance/sum(variance)*100, pct_cum = cumsum(pct))
  scores <- as_tibble(pca$x, rownames = "sample")
  colnames(coldata) <- c("sample", "STAG2_gt")
  input_pca <- merge(coldata, scores, by = "sample")
  
  p <- ggplot(input_pca, aes(PC1, PC2, colour = factor(STAG2_gt))) +
    geom_point() +
    labs(x = paste0("PC1 (", round(eig$pct[1], 1), "%)"),
         y = paste0("PC2 (", round(eig$pct[2], 1), "%)"))
  
  list(top_genes = top, eigenvalues = eig, merged = input_pca, plot = p)
}
res <- run_myeloid_pca(tpm_only_myeloid, All_myeloid_celllines_cleaned_with_Gt)
res$plot
