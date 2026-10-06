suppressMessages(library(edgeR))
set.seed(1)

run_da <- function(counts, samples, label) {
  samples$species <- factor(samples$species, levels = c("Pman", "Ppol"))
  samples$sex <- factor(samples$sex, levels = c("F", "M"))
  samples$rep <- factor(samples$rep)
  y <- DGEList(counts = counts, samples = samples)
  y$samples$lib.size <- colSums(counts)
  y <- calcNormFactors(y)
  design <- model.matrix(~ rep + species + sex, samples)
  y <- estimateDisp(y, design, trend = "none")
  fit <- glmQLFit(y, design)
  out <- list()
  for (term in c("speciesPpol", "sexM")) {
    res <- glmQLFTest(fit, coef = term)
    tt <- topTags(res, n = Inf, adjust.method = "BH")$table
    tt$cluster <- rownames(tt); tt$term <- term; tt$analysis <- label
    tt$fdr_holm <- p.adjust(tt$PValue, "holm")
    out[[term]] <- tt
  }
  do.call(rbind, out)
}

ct <- as.matrix(read.csv("results/cluster_by_animal_counts.csv", row.names = 1, check.names = FALSE))
meta <- read.csv("results/animal_metadata.csv", row.names = 1)
meta$species <- sub("P\\.man", "Pman", sub("P\\.pol", "Ppol", meta$species))
ct <- ct[, rownames(meta)]
resA <- run_da(ct, meta, "my_clusters")

md <- read.csv("data/neurons_clustered_metadata.csv")
md$animal <- paste(md$species, md$sex, md$rep, sep = "_")
ctB <- table(md$cluster, md$animal)
metaB <- unique(md[, c("animal", "species", "sex", "rep")]); rownames(metaB) <- metaB$animal
metaB$species <- ifelse(metaB$species == "BW", "Pman", "Ppol")
ctB <- unclass(ctB)[, rownames(metaB)]
resB <- run_da(ctB, metaB, "authors_clusters")

write.csv(rbind(resA, resB), "results/differential_abundance_edger.csv", row.names = FALSE)
for (a in c("authors_clusters", "my_clusters")) for (t in c("speciesPpol", "sexM")) {
  r <- rbind(resA, resB); r <- r[r$analysis == a & r$term == t, ]
  cat(sprintf("%s %s: FDR(BH)<0.05 = %d, FDR(Holm)<0.05 = %d, nominal p<0.05 = %d\n",
              a, t, sum(r$FDR < 0.05), sum(r$fdr_holm < 0.05), sum(r$PValue < 0.05)))
}
