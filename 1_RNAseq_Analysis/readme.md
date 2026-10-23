# DRPPM-EASY: Gene Expression Analysis in R (Without Shiny)

> **Hands-on tutorial:** Reproduce a simplified DRPPM-EASY-inspired gene expression workflow in an ordinary R session. This tutorial is an educational standalone implementation, **not** a claim that the code below is extracted verbatim from the Shiny application.

**Related project:** [DRPPM-EASY-Integration](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration)

## Learning objectives

By the end of this tutorial, you will be able to:

1. Import a gene-by-sample expression matrix and sample metadata.
2. Filter and log-transform normalized expression measurements.
3. Inspect sample relationships with principal component analysis (PCA).
4. Compare two experimental groups with `limma`.
5. Interpret differential expression statistics and create a volcano plot.
6. Visualize top-ranked genes in a heatmap.
7. Run a demonstration gene-set test using `limma::camera`.

**Estimated time:** 15–30 minutes (after package installation).  
**Requirements:** R, an internet connection for first-time package installation, and no Shiny server.

## 1. Download the tutorial files

Keep these three files in the same directory:

```text
EASY_GitHub_Tutorial/
├── README.md
├── EASY_gene_expression_tutorial.R
├── example_expression.tsv
├── example_metadata.tsv
└── results/                  # created when the script runs
```

The example data are **synthetic**, comprising 200 genes and 8 samples (4 Control and 4 Tumor). They are intended to demonstrate the analysis workflow, **not** to support biological conclusions.

If the tutorial files are added to the upstream GitHub repository, clone the repository and navigate to the folder where the files were placed. Otherwise, download this tutorial bundle and unzip it locally.

## 2. Install the R packages

Run once in the R console:

```r
install.packages(c("ggplot2", "pheatmap"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("limma")
```

Load packages to confirm they are installed:

```r
library(limma)
library(ggplot2)
library(pheatmap)
```

## 3. Understand the input files

### Expression matrix: `example_expression.tsv`

Rows are genes; columns are samples. The first column must be called `Gene`. The numeric values are **nonnegative, normalized, linear-scale expression measurements**, not raw read counts or already-log-transformed values.

| Gene | Control_1 | Control_2 | ... | Tumor_4 |
|---|---:|---:|---:|---:|
| Gene_001 | numeric | numeric | ... | numeric |
| Gene_002 | numeric | numeric | ... | numeric |

### Sample metadata: `example_metadata.tsv`

| Sample | Group |
|---|---|
| Control_1 | Control |
| Control_2 | Control |
| ... | ... |
| Tumor_4 | Tumor |

**Important:** Sample names must match the expression matrix columns. This tutorial uses exactly the group labels `Control` and `Tumor`, with at least two samples per group. Gene identifiers must be unique.

Preview both files in R:

```r
head(read.delim("example_expression.tsv", check.names = FALSE)[, 1:5])
read.delim("example_metadata.tsv")
```

## 4. Import and validate the data

```r
raw <- read.delim("example_expression.tsv", check.names = FALSE)
stopifnot(names(raw)[1] == "Gene", !anyDuplicated(raw$Gene))

expr <- as.matrix(raw[, -1, drop = FALSE])
rownames(expr) <- raw$Gene
storage.mode(expr) <- "numeric"

meta <- read.delim("example_metadata.tsv", stringsAsFactors = FALSE)
stopifnot(all(c("Sample", "Group") %in% names(meta)))
stopifnot(!anyDuplicated(meta$Sample))
stopifnot(all(meta$Sample %in% colnames(expr)))

# Align expression columns to metadata row order
expr <- expr[, meta$Sample, drop = FALSE]
stopifnot(all(is.finite(expr)), all(expr >= 0))

meta$Group <- factor(meta$Group, levels = c("Control", "Tumor"))
stopifnot(!anyNA(meta$Group), all(table(meta$Group) >= 2))

dim(expr)
table(meta$Group)
```

**Interpretation:** `dim(expr)` should report 200 rows and 8 columns for the supplied example. `table(meta$Group)` should show 4 samples in each group.

## 5. Filter and log-transform normalized expression

```r
# Retain genes with normalized expression > 1 in at least 3 samples
keep <- rowSums(expr > 1) >= 3
logexpr <- log2(expr[keep, , drop = FALSE] + 1)

dim(logexpr)
summary(as.vector(logexpr))
```

The `log2(x + 1)` transformation compresses high values and avoids taking the logarithm of zero. **Do not apply it again if your input is already log-transformed.** The filtering threshold is illustrative and should be tailored to the normalization method and study design.

## 6. Perform PCA for exploratory quality control

```r
pca <- prcomp(t(logexpr), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  Sample = meta$Sample, Group = meta$Group
)

p <- ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  theme_bw() +
  ggtitle("Sample PCA")

print(p)
ggsave("results/PCA.png", p, width = 7, height = 5, dpi = 160)
```

**How to read the PCA:** Samples close together have similar global expression profiles along the displayed components. Look for outliers, batch effects, and separation between groups. Group separation is exploratory, **not** evidence by itself that individual genes are differentially expressed.

## 7. Differential expression with `limma`

Here, the contrast is **Tumor minus Control**. A positive `logFC` indicates higher expression in Tumor; a negative value indicates lower expression.

```r
design <- model.matrix(~ 0 + Group, data = meta)
colnames(design) <- levels(meta$Group)

fit <- lmFit(logexpr, design)
contrast <- makeContrasts(Tumor - Control, levels = design)
fit <- eBayes(contrasts.fit(fit, contrast))

deg <- topTable(fit, number = Inf, sort.by = "P", adjust.method = "BH")
deg$Gene <- rownames(deg)
deg <- deg[, c("Gene", setdiff(names(deg), "Gene"))]

deg$Significance <- ifelse(
  deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 1,
  ifelse(deg$logFC > 0, "Up", "Down"),
  "Not significant"
)

write.csv(deg, "results/differential_expression.csv", row.names = FALSE)
head(deg)
table(deg$Significance)
```

| Column | Meaning |
|---|---|
| `Gene` | Gene identifier |
| `logFC` | Estimated Tumor–Control difference on the log2-expression scale |
| `AveExpr` | Average log2 expression |
| `t` | Moderated t-statistic |
| `P.Value` | Unadjusted P-value |
| `adj.P.Val` | Benjamini–Hochberg adjusted P-value (FDR) |
| `B` | Log-odds of differential expression under the fitted model |
| `Significance` | Tutorial label using FDR < 0.05 and \|logFC\| ≥ 1 |

**Statistical note:** This simple design does not account for pairing, batches, or covariates. In real studies, include appropriate variables in the design and use methods suited to the measurement scale. The fold-change cutoff here is a *post hoc classification threshold*, not a formal minimum-effect hypothesis test.

## 8. Create a volcano plot

```r
v <- ggplot(
  deg,
  aes(logFC, -log10(pmax(adj.P.Val, 1e-300)), color = Significance)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Tumor vs Control",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  )

print(v)
ggsave("results/volcano.png", v, width = 7, height = 5, dpi = 160)
```

**Interpretation:** Points to the right have higher expression in Tumor; points to the left have lower expression. Points higher on the plot have smaller adjusted P-values. This example uses adjusted P-values on the y-axis (some volcano plots use raw P-values).

## 9. Plot a heatmap of the top 30 genes

```r
top <- head(deg$Gene, 30)
ann <- data.frame(Group = meta$Group, row.names = meta$Sample)

pheatmap(
  logexpr[top, , drop = FALSE],
  scale = "row",
  annotation_col = ann,
  filename = "results/heatmap.png",
  width = 8, height = 7,
  main = "Top differentially expressed genes"
)
```

**Interpretation:** Row scaling emphasizes relative variation across samples *within each gene*. Colors should not be interpreted as directly comparable absolute expression between different genes. The heatmap uses the top-ranked genes by P-value; they are not necessarily all significant.

## 10. Demonstrate gene-set enrichment with `camera`

The tutorial uses **synthetic demonstration gene sets** so that the code works offline after package installation. These sets are **not** curated pathways. For publication-quality biology, substitute gene sets from MSigDB, GO, Reactome, or another validated source, and ensure identifiers match.

```r
sets <- list(
  Demo_Upregulated = rownames(logexpr)[1:25],
  Demo_Downregulated = rownames(logexpr)[26:50],
  Demo_Background = rownames(logexpr)[51:90]
)

indices <- lapply(sets, function(genes) {
  which(rownames(logexpr) %in% genes)
})

enrich <- camera(
  logexpr,
  index = indices,
  design = design,
  contrast = contrast
)
enrich$GeneSet <- rownames(enrich)
write.csv(enrich, "results/demo_gene_set_enrichment.csv", row.names = FALSE)
print(enrich)
```

**Interpretation:** `camera` performs a competitive gene-set test accounting for inter-gene correlation. It is **not identical to classic ranked GSEA or ssGSEA**. Examine direction, gene-set size, and FDR before drawing biological conclusions.

## 11. Run the complete workflow in one command

In RStudio, set the working directory to the folder containing the tutorial files, then execute:

```r
setwd("/path/to/EASY_GitHub_Tutorial")  # replace with your local path
source("EASY_gene_expression_tutorial.R")
```

The script creates `results/` automatically. It saves:

| Output | Description |
|---|---|
| `results/PCA.png` | PCA of the expression matrix |
| `results/differential_expression.csv` | Gene-level differential expression statistics |
| `results/volcano.png` | Differential expression volcano plot |
| `results/heatmap.png` | Heatmap of the top 30 ranked genes |
| `results/demo_gene_set_enrichment.csv` | Demonstration gene-set results |

## 12. Apply the tutorial to your own data

Replace `example_expression.tsv` and `example_metadata.tsv` with files in the same structure, or edit the two filenames in the script. Confirm that:

- Genes are unique row identifiers, and sample names match exactly.
- Expression values are normalized, nonnegative, numeric, and on a **linear** scale.
- Groups are labeled `Control` and `Tumor`, or update the factor levels and contrast in the script.
- Each group has biological replicates; include relevant covariates for your experimental design.
- Genes and pathway sets use compatible identifier systems.

### What if my input is raw RNA-seq counts?

**Do not** run the simple `log2(counts + 1)` + `limma` analysis above as a substitute for count-aware normalization and precision modeling. For raw integer counts, use `edgeR` filtering and TMM normalization followed by `voom` + `limma`, or use `DESeq2`. The supplied standalone script is designed for normalized expression data.

## Troubleshooting

| Problem | Likely cause / solution |
|---|---|
| `cannot open file` | Set the working directory to the folder containing the `.tsv` files. |
| `Install package: limma` | Install `BiocManager`, then run `BiocManager::install("limma")`. |
| `all(meta$Sample %in% colnames(expr)) is not TRUE` | Sample IDs differ between metadata and expression columns. |
| `all(is.finite(expr)) is not TRUE` | Missing, infinite, or nonnumeric expression values need cleaning. |
| `all(expr >= 0) is not TRUE` | Input may already be centered/scaled or otherwise unsuitable for the log transform. |
| No significant genes | Check group labels, sample size, variability, normalization, and design. This can be a valid result. |
| PCA or heatmap looks unexpected | Inspect outliers, sample swaps, batch effects, and whether input is already log-transformed. |

## Further reading

- [DRPPM-EASY-Integration GitHub repository](https://github.com/shawlab-moffitt/DRPPM-EASY-Integration)
- [limma Bioconductor documentation](https://bioconductor.org/packages/limma/)
- [edgeR Bioconductor documentation](https://bioconductor.org/packages/edgeR/)
- [MSigDB gene-set collections](https://www.gsea-msigdb.org/gsea/msigdb/)

---

**Scope and limitations:** This walkthrough illustrates PCA, two-group differential expression, visualization, and a competitive gene-set test using simulated data. It does not implement every feature of the DRPPM-EASY Shiny application, nor does it validate the results against the repository's internal analytical functions.
