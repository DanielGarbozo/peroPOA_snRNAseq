suppressMessages({library(edgeR); library(Matrix)})
set.seed(1)
MODE <- Sys.getenv("SEXDE_MODE", "strict")
MIN_NUCLEI <- if (MODE == "strict") 10 else 1
FDR_CUT <- 0.05

run_one <- function(tag, label) {
  pb <- readMM(sprintf("data/pseudobulk_counts%s.mtx", tag))
  genes <- readLines("data/pseudobulk_genes.csv")
  meta <- read.csv(sprintf("data/pseudobulk_meta%s.csv", tag))
  rownames(pb) <- genes
  res <- list()
  for (sp in c("P.man", "P.pol")) {
    for (cl in unique(meta$cluster)) {
      idx <- which(meta$species == sp & meta$cluster == cl & meta$n_nuclei >= MIN_NUCLEI)
      m <- meta[idx, ]
      if (length(unique(m$sex)) < 2 || sum(m$sex == "F") < 2 || sum(m$sex == "M") < 2) next
      m$rep <- factor(m$rep); m$sex <- factor(m$sex, levels = c("F", "M"))
      design <- model.matrix(~ rep + sex, m)
      if (qr(design)$rank < ncol(design) || nrow(m) - ncol(design) < 2) next
      y <- DGEList(counts = as.matrix(pb[, idx]))
      keep <- if (MODE == "strict") filterByExpr(y, design) else rowSums(y$counts) > 0
      y <- y[keep, , keep.lib.sizes = FALSE]
      y <- calcNormFactors(y)
      y <- estimateDisp(y, design)
      fit <- glmQLFit(y, design)
      tt <- topTags(glmQLFTest(fit, coef = "sexM"), n = Inf, adjust.method = "BH")$table
      tt$gene <- rownames(tt); tt$cluster <- cl; tt$species <- sp; tt$n_samples <- nrow(m)
      res[[paste(sp, cl)]] <- tt
    }
  }
  out <- do.call(rbind, res); out$analysis <- label
  out
}

TAGS <- Sys.getenv("SEXDE_TAGS", "")
if (nzchar(TAGS)) {
  all <- do.call(rbind, lapply(strsplit(TAGS, ";")[[1]], function(t) run_one(t, sub("^_", "", t))))
} else {
  all <- rbind(run_one("", "my_clusters"), run_one("_auth", "authors_clusters"))
}
sfx <- if (MODE == "strict") "" else paste0("_", MODE)
if (!nzchar(TAGS)) {
  write.csv(all, sprintf("data/sex_de_edger%s.csv", sfx), row.names = FALSE)
  write.csv(all[all$FDR < FDR_CUT, ], sprintf("results/sex_de_edger%s_sig.csv", sfx), row.names = FALSE)
}

summ <- list()
for (an in unique(all$analysis)) for (sp in c("P.man", "P.pol")) {
  d <- all[all$analysis == an & all$species == sp, ]
  sig <- d[d$FDR < FDR_CUT, ]
  top <- do.call(rbind, lapply(split(sig, sig$gene), function(x) x[which.max(x$logCPM), ]))
  summ[[paste(an, sp)]] <- data.frame(analysis = an, species = sp,
    clusters_tested = length(unique(d$cluster)),
    sig_gene_cluster_pairs = nrow(sig),
    genes_male_biased = sum(top$logFC > 0), genes_female_biased = sum(top$logFC < 0))
}
summ <- do.call(rbind, summ); rownames(summ) <- NULL
summ$mode <- MODE
out_file <- if (nzchar(TAGS)) sprintf("results/robustness/edger_summary%s.csv", Sys.getenv("SEXDE_OUT", "_custom")) else sprintf("results/sex_de_summary%s.csv", sfx)
write.csv(summ, out_file, row.names = FALSE)
print(summ)
