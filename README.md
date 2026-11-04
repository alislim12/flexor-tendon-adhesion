# MMP Targets in Flexor Tendon Adhesion

Independent study (Bio 593) investigating matrix metalloproteinase (MMP)
family genes in flexor tendon adhesion and hand fibrosis, using public
transcriptomic data.

## Summary
TIMP3 (a master MMP inhibitor) was identified as a candidate anti-adhesion
gene and validated across four independent comparisons:
- Mouse flexor tendon, injured vs. uninjured (down)
- Mouse flexor tendon, PAI-1 knockout vs. wild-type (up when scarring reduced)
- Human Dupuytren's tissue, bulk (down)
- Human Dupuytren's fibroblasts, single-cell (down)

## Data (downloaded from GEO, not included in repo)
- GSE125026 — mouse flexor tendon RNA-seq (bulk)
- GSE75152 — human Dupuytren's microarray (bulk)
- GSE173252 — human Dupuytren's single-cell RNA-seq

## Scripts
- `analysis.R` — bulk RNA-seq analysis (DESeq2)
- `singlecell.R` — single-cell analysis (Seurat)

## Results
CSV files of differential-expression results and the volcano/comparison figures.