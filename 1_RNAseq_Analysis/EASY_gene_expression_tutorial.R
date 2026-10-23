# DRPPM-EASY-inspired gene expression tutorial (standalone R, no Shiny)
# Input: normalized nonnegative expression values, NOT raw RNA-seq counts.
# Run from the directory containing example_expression.tsv and example_metadata.tsv.
# Install once: install.packages(c('ggplot2','pheatmap'))
# if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')
# BiocManager::install('limma')

for (pkg in c('limma','ggplot2','pheatmap')) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop('Install package: ', pkg)
}
dir.create('results', showWarnings=FALSE)

# 1. Read and align expression and metadata
raw <- read.delim('example_expression.tsv', check.names=FALSE)
stopifnot(names(raw)[1] == 'Gene', !anyDuplicated(raw$Gene))
expr <- as.matrix(raw[, -1, drop=FALSE])
rownames(expr) <- raw$Gene
storage.mode(expr) <- 'numeric'
meta <- read.delim('example_metadata.tsv', stringsAsFactors=FALSE)
stopifnot(all(c('Sample','Group') %in% names(meta)),
          !anyDuplicated(meta$Sample), all(meta$Sample %in% colnames(expr)))
expr <- expr[, meta$Sample, drop=FALSE]
stopifnot(all(is.finite(expr)), all(expr >= 0))
meta$Group <- factor(meta$Group, levels=c('Control','Tumor'))
stopifnot(!anyNA(meta$Group), all(table(meta$Group) >= 2))

# 2. Log2 transform and filter low-expression genes
# Demo data are normalized linear-scale values; do NOT log2-transform twice.
keep <- rowSums(expr > 1) >= 3
logexpr <- log2(expr[keep, , drop=FALSE] + 1)

# 3. PCA quality control
pca <- prcomp(t(logexpr), center=TRUE, scale.=FALSE)
pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],
                     Sample=meta$Sample, Group=meta$Group)
p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color=Group, label=Sample)) +
  ggplot2::geom_point(size=3) + ggplot2::geom_text(vjust=-0.7, size=3) +
  ggplot2::theme_bw() + ggplot2::ggtitle('Sample PCA')
ggplot2::ggsave('results/PCA.png', p, width=7, height=5, dpi=160)

# 4. Differential expression using limma moderated t-tests
# Tumor - Control: positive logFC = higher in Tumor.
design <- model.matrix(~0 + Group, data=meta)
colnames(design) <- levels(meta$Group)
fit <- limma::lmFit(logexpr, design)
contrast <- limma::makeContrasts(Tumor-Control, levels=design)
fit <- limma::eBayes(limma::contrasts.fit(fit, contrast))
deg <- limma::topTable(fit, number=Inf, sort.by='P', adjust.method='BH')
deg$Gene <- rownames(deg)
deg <- deg[, c('Gene', setdiff(names(deg), 'Gene'))]
deg$Significance <- ifelse(deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 1,
                           ifelse(deg$logFC > 0, 'Up', 'Down'), 'Not significant')
write.csv(deg, 'results/differential_expression.csv', row.names=FALSE)
print(table(deg$Significance))

# 5. Volcano plot
v <- ggplot2::ggplot(deg, ggplot2::aes(logFC, -log10(pmax(adj.P.Val, 1e-300)),
                                          color=Significance)) +
  ggplot2::geom_point(alpha=0.7, size=1.5) +
  ggplot2::geom_vline(xintercept=c(-1,1), linetype='dashed') +
  ggplot2::geom_hline(yintercept=-log10(0.05), linetype='dashed') +
  ggplot2::theme_bw() + ggplot2::labs(title='Tumor vs Control', x='log2 fold change',
                                      y='-log10 adjusted p-value')
ggplot2::ggsave('results/volcano.png', v, width=7, height=5, dpi=160)

# 6. Heatmap of 30 genes with smallest adjusted p-values
n <- min(30, nrow(deg))
top <- deg$Gene[seq_len(n)]
ann <- data.frame(Group=meta$Group, row.names=meta$Sample)
pheatmap::pheatmap(logexpr[top, , drop=FALSE], scale='row',
                   annotation_col=ann, filename='results/heatmap.png',
                   width=8, height=7, main='Top differentially expressed genes')

# 7. Illustrative gene-set enrichment with limma::camera
# These are SYNTHETIC example sets, NOT validated biological pathways.
# Replace with curated MSigDB/GO gene sets for real analyses.
sets <- list(Demo_Upregulated=rownames(logexpr)[1:25],
             Demo_Downregulated=rownames(logexpr)[26:50],
             Demo_Background=rownames(logexpr)[51:90])
indices <- lapply(sets, function(x) which(rownames(logexpr) %in% x))
enrich <- limma::camera(logexpr, index=indices, design=design, contrast=contrast)
enrich$GeneSet <- rownames(enrich)
write.csv(enrich, 'results/demo_gene_set_enrichment.csv', row.names=FALSE)
print(enrich)
cat('\nSaved outputs to results/\n')

# For raw RNA-seq counts instead: use edgeR TMM normalization + voom/limma
# (or DESeq2), rather than applying this normalized-expression workflow directly.
