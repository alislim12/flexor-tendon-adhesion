# ── Single-cell analysis: GSE173252 Dupuytren (TIMP3 in fibroblasts) ──
library(Seurat)
library(R.utils)

# ---- Load the pre-processed Seurat object ----
# First time only: the GEO file was double-gzipped, so it was unpacked twice.
# After the first run, just load the saved copy below instead.
# R.utils::gunzip("GSE173252/GSE173252_dd.rds.gz", destname="GSE173252/dd_unzipped",
#                 remove=FALSE, overwrite=TRUE)
# R.utils::gunzip("GSE173252/dd_unzipped", destname="GSE173252/dd_final.rds",
#                 remove=FALSE, overwrite=TRUE)
# dd <- readRDS("GSE173252/dd_final.rds")
# dd <- UpdateSeuratObject(dd)
# saveRDS(dd, "GSE173252/dd_seurat.rds")

dd <- readRDS("GSE173252/dd_seurat.rds")   # the clean saved copy

# ---- TIMP3 in Dupuytren vs normal-fascia fibroblasts ----
fibro <- subset(dd, subset = annotation == "Mesenchyme")
fibro <- subset(fibro, subset = source %in% c("DD", "SF"))

DefaultAssay(fibro) <- "RNA"
fibro <- NormalizeData(fibro)
Idents(fibro) <- "source"

de <- FindMarkers(fibro, ident.1 = "DD", ident.2 = "SF",
                  features = c("TIMP3", "MMP12", "MMP2", "MMP14"))
print(de)

write.csv(de, "GSE173252_fibroblast_DD_vs_SF_TIMP3.csv")







# ── Single-cell Step 3: Full DE across all genes (DD vs SF fibroblasts) ──
# No 'features' argument this time = test every gene
de_all <- FindMarkers(fibro, ident.1 = "DD", ident.2 = "SF",
                      min.pct = 0.1, logfc.threshold = 0.25)

cat("Genes tested / passing thresholds:", nrow(de_all), "\n\n")
cat("Top 15 by significance:\n")
print(head(de_all[order(de_all$p_val_adj), ], 15))

# Where does TIMP3 rank in the full list?
cat("\nTIMP3 in the full results:\n")
print(de_all["TIMP3", ])

write.csv(de_all, "GSE173252_fibroblast_DD_vs_SF_allgenes.csv")











# ── Single-cell Step 4: Volcano plot of DD vs SF fibroblasts ──────
de_all$gene <- rownames(de_all)

# Cap p-values that are exactly 0 so they can be plotted (log of 0 = infinity)
de_all$padj_plot <- de_all$p_val_adj
de_all$padj_plot[de_all$padj_plot == 0] <- 1e-300
de_all$neglog10 <- -log10(de_all$padj_plot)

# Base plot: all genes in grey
plot(de_all$avg_log2FC, de_all$neglog10,
     pch = 20, col = "grey70", cex = 0.5,
     xlab = "log2 fold-change (DD vs normal fascia)",
     ylab = "-log10 adjusted p-value",
     main = "Dupuytren vs normal fibroblasts (single-cell)")
abline(v = c(-1, 1), lty = 2, col = "grey60")

# Highlight genes of interest in red + label them
genes_label <- c("TIMP3", "MMP2", "MMP14", "POSTN", "ADAM12", "COL5A1", "TNC")
hits <- de_all$gene %in% genes_label
points(de_all$avg_log2FC[hits], de_all$neglog10[hits],
       col = "red", pch = 20, cex = 1.1)
text(de_all$avg_log2FC[hits], de_all$neglog10[hits],
     labels = de_all$gene[hits], pos = 3, cex = 0.8, col = "red")










# ── Single-cell Step 5: Pathway enrichment (Enrichr) ──────────────
library(enrichR)

# Split changed genes by direction (using a modest fold-change cutoff)
up_genes   <- rownames(de_all)[de_all$avg_log2FC >  0.5]
down_genes <- rownames(de_all)[de_all$avg_log2FC < -0.5]
cat("Up in DD:", length(up_genes), " | Down in DD:", length(down_genes), "\n")

# Choose a few well-known gene-set libraries
dbs_use <- c("GO_Biological_Process_2023",
             "Reactome_2022",
             "MSigDB_Hallmark_2020")

# Run enrichment for each direction
enr_up   <- enrichr(up_genes,   dbs_use)
enr_down <- enrichr(down_genes, dbs_use)

# Show the top pathways for UP genes (the fibrosis program)
cat("\n=== TOP PATHWAYS: UP in Dupuytren ===\n")
print(head(enr_up[["GO_Biological_Process_2023"]][, c("Term","Overlap","Adjusted.P.value")], 10))

# And for DOWN genes (where TIMP3 lives)
cat("\n=== TOP PATHWAYS: DOWN in Dupuytren ===\n")
print(head(enr_down[["GO_Biological_Process_2023"]][, c("Term","Overlap","Adjusted.P.value")], 10))