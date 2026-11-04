# Flexor Tendon Adhesion — MMP target analysis (GSE125026)
library(GEOquery)
library(DESeq2)

# Load raw counts (genes x samples)
counts <- read.delim("GSE125026/GSE125026_deSeq2_counts.txt.gz",
                     row.names = 1, check.names = FALSE)
# Download sample metadata and build 'pheno' (needed for sample labels)
gse   <- getGEO("GSE125026", GSEMatrix = TRUE)
gse   <- gse[[1]]
pheno <- pData(gse)

# ── Step 3: Build sample sheet from the TRUE titles, not column names ──
# Parse genotype + condition out of the reliable 'title' field
titles <- pheno$title

coldata <- data.frame(
  row.names = colnames(counts),          # keep count-column order
  title     = titles,
  genotype  = ifelse(grepl("^WT",  titles), "WT", "PAI1"),
  condition = dplyr::case_when(
    grepl("UNINJ", titles) ~ "uninjured",
    grepl("2wPL",  titles) ~ "inj_2wk",
    grepl("4wPL",  titles) ~ "inj_4wk"
  )
)

# Sanity check: how many samples fall in each genotype x condition cell?
print(table(coldata$genotype, coldata$condition))

# ── Step 4: Differential expression — WT injured (2wk) vs WT uninjured ──

# 1. Keep only the WT 2-week and WT uninjured samples (6 total)
keep <- coldata$genotype == "WT" &
  coldata$condition %in% c("inj_2wk", "uninjured")

counts_sub  <- counts[, keep]
coldata_sub <- coldata[keep, ]

# Make 'uninjured' the reference level, so results read as injured-vs-uninjured
coldata_sub$condition <- factor(coldata_sub$condition,
                                levels = c("uninjured", "inj_2wk"))

cat("Samples going into this comparison:\n")
print(coldata_sub[, c("title", "genotype", "condition")])

# 2. Build the DESeq2 object and run the analysis
dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                              colData   = coldata_sub,
                              design    = ~ condition)

# Drop genes with almost no reads (speeds it up, removes noise)
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- DESeq(dds)       # this is the actual statistical run

# 3. Extract results
res <- results(dds, contrast = c("condition", "inj_2wk", "uninjured"))
summary(res)




# ── Step 5: Focus on the MMP family and their regulators ──────────
res_ordered <- res[order(res$padj), ]   # sort results by significance
res_df <- as.data.frame(res_ordered)
res_df$gene <- rownames(res_df)

# Mouse gene symbols: MMPs are "Mmp...", their inhibitors are "Timp..."
# (Adamts = related matrix-remodeling proteases; included as context)
mmp_pattern <- "^Mmp|^Timp|^Adamts"

mmp_hits <- res_df[grepl(mmp_pattern, res_df$gene), ]

# Show them ranked by significance, with the columns that matter
mmp_table <- mmp_hits[order(mmp_hits$padj),
                      c("gene", "log2FoldChange", "pvalue", "padj", "baseMean")]

cat("MMP-family genes detected:", nrow(mmp_table), "\n\n")
print(mmp_table)






# ── Step 6: PAI-1 KO vs WT, both injured at 2 weeks ───────────────
# Question: what changes when the scar-related gene (PAI-1) is disrupted?

# 1. Keep only the 2-week injured samples from BOTH genotypes (6 total)
keep2 <- coldata$condition == "inj_2wk"

counts_sub2  <- counts[, keep2]
coldata_sub2 <- coldata[keep2, ]

# WT is the reference, so results read as "PAI1 relative to WT"
coldata_sub2$genotype <- factor(coldata_sub2$genotype,
                                levels = c("WT", "PAI1"))

cat("Samples in this comparison:\n")
print(coldata_sub2[, c("title", "genotype", "condition")])

# 2. Run DESeq2 (design is now ~ genotype)
dds2 <- DESeqDataSetFromMatrix(countData = counts_sub2,
                               colData   = coldata_sub2,
                               design    = ~ genotype)
dds2 <- dds2[rowSums(counts(dds2)) >= 10, ]
dds2 <- DESeq(dds2)

res2 <- results(dds2, contrast = c("genotype", "PAI1", "WT"))
summary(res2)


















# ── Step 6b: MMP family in the PAI-1 KO vs WT comparison ──────────
res2_df <- as.data.frame(res2)
res2_df$gene <- rownames(res2_df)

mmp_hits2 <- res2_df[grepl("^Mmp|^Timp|^Adamts", res2_df$gene), ]
mmp_table2 <- mmp_hits2[order(mmp_hits2$padj),
                        c("gene", "log2FoldChange", "pvalue", "padj", "baseMean")]

cat("MMP-family genes (PAI-1 KO vs WT, 2wk):", nrow(mmp_table2), "\n\n")
print(mmp_table2)

# Save it alongside the first one
write.csv(mmp_table2, "MMP_results_PAI1_vs_WT_2wk.csv", row.names = FALSE)













# ── Step 7: Cross-reference the two MMP comparisons ───────────────
# Merge the injury signature with the knockout effect, by gene
merged <- merge(
  mmp_table[,  c("gene", "log2FoldChange", "padj")],
  mmp_table2[, c("gene", "log2FoldChange", "padj")],
  by = "gene",
  suffixes = c("_injury", "_KO")
)

# Candidates: significant in BOTH comparisons
both_sig <- subset(merged, padj_injury < 0.1 & padj_KO < 0.1)

# Rank by combined strength (sum of -log10 padj)
both_sig$score <- -log10(both_sig$padj_injury) + -log10(both_sig$padj_KO)
both_sig <- both_sig[order(-both_sig$score), ]

cat("MMP genes significant in BOTH comparisons:", nrow(both_sig), "\n\n")
print(both_sig)

write.csv(both_sig, "MMP_candidates_overlap.csv", row.names = FALSE)















# ── Step 8: Visualize — injury effect vs knockout effect ──────────
plot_df <- merge(
  mmp_table[,  c("gene", "log2FoldChange")],
  mmp_table2[, c("gene", "log2FoldChange")],
  by = "gene", suffixes = c("_injury", "_KO")
)

plot(plot_df$log2FoldChange_injury, plot_df$log2FoldChange_KO,
     xlab = "log2FC: injury vs uninjured (WT)",
     ylab = "log2FC: PAI-1 KO vs WT (injured)",
     main = "MMP-family genes: scarring vs anti-scarring",
     pch = 19, col = "grey50")
abline(h = 0, v = 0, lty = 2, col = "grey70")

# Label the two overlap candidates in red
hits <- plot_df$gene %in% c("Timp3", "Mmp12")
points(plot_df$log2FoldChange_injury[hits],
       plot_df$log2FoldChange_KO[hits], col = "red", pch = 19)
text(plot_df$log2FoldChange_injury[hits],
     plot_df$log2FoldChange_KO[hits],
     labels = plot_df$gene[hits], pos = 4, col = "red")















# ── Step 9: Human cross-check — Dupuytren (GSE75152) ──────────────
# Do TIMP3 / MMP12 show up in human hand fibrosis too?

gse_human <- getGEO("GSE75152", GSEMatrix = TRUE)
gse_human <- gse_human[[1]]

cat("Samples:", ncol(gse_human), "\n")
cat("Probes/genes:", nrow(gse_human), "\n\n")

# What do the sample groups look like? (disease vs control)
print(pData(gse_human)[, c("title", "geo_accession")])








# ── Step 9b: Test TIMP3 / MMP12 in human Dupuytren (GSE75152) ─────
expr  <- exprs(gse_human)                       # expression matrix
fdat  <- fData(gse_human)                        # probe -> gene annotation
group <- ifelse(grepl("patient", pData(gse_human)$title),
                "Dupuytren", "control")
group <- factor(group, levels = c("control", "Dupuytren"))

# Helper: for a gene symbol, average its probes per sample, then
# compare Dupuytren vs control with a t-test
test_gene <- function(sym) {
  probes <- rownames(fdat)[fdat$Symbol == sym]
  probes <- probes[probes %in% rownames(expr)]
  if (length(probes) == 0) return(cat(sym, ": not found on this array\n"))
  
  vals <- colMeans(expr[probes, , drop = FALSE])   # avg across probes
  tt   <- t.test(vals ~ group)
  log2fc <- mean(vals[group=="Dupuytren"]) - mean(vals[group=="control"])
  
  cat(sprintf("%-6s | probes: %d | log2FC (Dupuytren vs control): %+.2f | p = %.3g\n",
              sym, length(probes), log2fc, tt$p.value))
}

cat("Human Dupuytren vs control:\n")
test_gene("TIMP3")
test_gene("MMP12")