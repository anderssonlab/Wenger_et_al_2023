library(limma)
library(DESeq2) 
library(tximport)
library(clusterProfiler)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(BiocParallel)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(pheatmap)
library(openxlsx)
library(Rsubread)
library(ggthemr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggpubr)
library(ggsci)
library(edgeR)
library(RColorBrewer)
library(latex2exp)
library(clusterProfiler)
library(ReactomePA)
library(reactome.db)
library(marray)
library(deco)
library(parallel)
library(reshape2)
library(sva)
library(Rsubread)
library(ggvenn)
library(ecoflux)
library(sva)
library(seqsetvis)

wd.dir <- "/isdata/alab/people/nalcaraz/Projects/KU/rep_chromatin_TC/publication/Wenger_et_al_2022"
setwd(wd.dir)
source("utils/utils.R")

plots.dir <- "plots"

data.dir <- "data"
rnaseq.data.dir <- file.path(data.dir, "RNA_seq")
external.data.dir <- file.path(data.dir, "external")

annot.version <- "gencode.vM23"
assembly <- "mm10"

chrom.sizes.file <- file.path(data.dir, paste0(assembly, ".chrom.sizes"))
chrom.sizes.df <- read.csv2(file = chrom.sizes.file, header = FALSE, col.names = c("chrom", "length"), sep = "\t")
chroms <- c(paste0("chr",1:19),"chrX","chrY")

annot.data.gtf <- file.path(external.data.dir,paste0(annot.version,".annotation.gtf.gz"))
annot.sqlite.file <- gsub(".gtf.gz",".sqlite", annot.data.gtf)

if (!file.exists(annot.sqlite.file)) {
  gencode.txdb <- makeTxDbFromGFF(file = annot.data.gtf, format="gtf", chrominfo=chrom.sizes.df, 
                                  dataSource=paste("arke:", annot.data.gtf, sep = ""),
                                  organism="Mus musculus")
  saveDb(gencode.txdb, file = annot.sqlite.file)
} else {
  gencode.txdb <- loadDb(annot.sqlite.file)
}

load(file.path(external.data.dir, "Processed_external.RData"))

gene.promoters.gr <- promoters(gencode.txdb, upstream=1000, 
                               downstream=1000, use.names=TRUE)

gene.promoters.gr  <- subset(gene.promoters.gr, seqnames %in% chroms)

k <- keys(gencode.txdb, keytype = "GENEID") 
df <- AnnotationDbi::select(gencode.txdb, keys = k, columns = "TXNAME", keytype = "GENEID")
tx2gene <- df[,2:1]
tx2gene$GENEID <- strsplit2(tx2gene$GENEID, "\\.")[,1]
rownames(tx2gene) <- tx2gene$TXNAME

gene.promoters.gr$gene.id <- tx2gene[gene.promoters.gr$tx_name,]$GENEID
gene.promoters.gr$CpG <- overlapsAny(gene.promoters.gr, 
                                     cpg.islands.gr)

gene.cpg.df <- gene.promoters.gr %>% as.data.frame() %>%
  group_by(gene.id, CpG) %>% tally() %>%
  mutate(CpG.frac = n / sum(n)) %>% 
  filter(CpG) %>% column_to_rownames(var = "gene.id") %>%
  as.data.frame()

genes.gr <- genes(gencode.txdb)
genes.gr <- subset(genes.gr, seqnames %in% chroms)

nrst.gene <- nearest(genes.gr, RT.gr)
names(genes.gr) <- substr(names(genes.gr),1,18)

genes.gr$RT <- RT.gr$RT.class[nrst.gene]

num.cores <- min(24,detectCores())

workers <- register(MulticoreParam(workers = num.cores))
register(MulticoreParam(workers = workers))


### ================= MCM2-2A Gene expression analyses =========================

gene.annot.df <- read.csv(file.path(rnaseq.data.dir, 
                                    "gene_and_TE_annotations.tsv.gz"),
                          header = TRUE,
                          sep = "\t")

clone.labels <- c("410" = "WT #1","506" = "WT #2","508" = "WT #3","509" = "WT #4",
                  "422" = "WT #5","503" = "WT #6","507" = "WT #7","510" = "WT #8",
                  "439" = "MCM2-2A #1", "442" = "MCM2-2A #2", "438" = "MCM2-2A #3", "421" = "MCM2-2A #4",
                  "443" = "MCM2-2A #5", "417" = "MCM2-2A #6", "418" = "MCM2-2A #7", "419" = "MCM2-2A #8")

rownames(gene.annot.df) <- gene.annot.df$ID
samples.annot.df <- read.csv(file.path(rnaseq.data.dir, "MCM2_2A_samples.tsv"), 
                             header = TRUE, sep = "\t")

rownames(samples.annot.df) <- samples.annot.df$ID

samples.annot.df$clone <- factor(samples.annot.df$clone, 
                                levels = c("410","422","503","506","507","508","509","510",
                                           "417","418","419","421","438","439","442","443"))

samples.annot.df$condition <- factor(samples.annot.df$condition, levels = c("WT","MCM2_2A"))

samples.annot.df$run <- factor(samples.annot.df$run, levels = c("AUG","JAN","MAR"))

samples.annot.df$labels <- clone.labels[as.character(samples.annot.df$clone)]


gene.counts.df <- read.csv(file.path(rnaseq.data.dir, "MCM2_gene_count_matrix.tsv.gz"),
                            header = TRUE, sep = "\t")

rownames(gene.counts.df) <- gene.counts.df$gene.id
gene.counts.df <- gene.counts.df[,-1]

log2FC.cutoff <- log2(1.5)
padj.cutoff <- 0.01

mcm2.dds <- DESeqDataSetFromMatrix(countData = as.matrix(gene.counts.df),
                                   colData = samples.annot.df[colnames(gene.counts.df),],
                                   rowData = gene.annot.df[rownames(gene.counts.df),],
                                   design = ~ run + condition)

mcm2.design <- model.matrix(~ run + clone, samples.annot.df)
to.keep <- filterByExpr(assay(mcm2.dds),design = mcm2.design, min.count = 1)

mcm2.ddq <- DESeq(mcm2.dds, parallel = TRUE)
mcm2.res <- results(mcm2.ddq, name = "condition_MCM2_2A_vs_WT")
mcm2.res.df <- as.data.frame(mcm2.res)
mcm2.res.df <- cbind(mcm2.res.df, as.data.frame(rowData(mcm2.dds)))

mcm2.gene.de.list <- list("Up" = filter(mcm2.res.df, 
                                        log2FoldChange >= log2FC.cutoff,
                                        padj < padj.cutoff)$ensembl_gene_id,
                          "Down" = filter(mcm2.res.df, 
                                          log2FoldChange <= -log2FC.cutoff,
                                          padj < padj.cutoff)$ensembl_gene_id)

mcm2.vst <- vst(assay(mcm2.dds,"counts"), blind = TRUE)

pca <- prcomp(t(mcm2.vst))

percentVar <- pca$sdev^2/sum(pca$sdev^2)
pca.df <- cbind(samples.annot.df[rownames(pca$x),], pca$x)

wt.clones <- unique(filter(samples.annot.df , condition == "WT")$clone)
wt.colors <- c(rev(brewer.pal(9, "Blues"))[1:4],rev(brewer.pal(9, "Greens"))[1:4])
names(wt.colors) <- as.character(wt.clones)
mut.clones <- unique(filter(samples.annot.df , condition == "MCM2_2A")$clone)
mut.colors <- c(rev(brewer.pal(9, "Oranges"))[1:4],rev(brewer.pal(9, "Reds"))[1:4])
names(mut.colors) <- as.character(mut.clones)
cond.colors <- c("MCM2_2A" = "#D55E00", "WT" = "#0072B2")


pca.x <- "PC1"
pca.y <- "PC2"

names(percentVar) <- colnames(pca$x)

pca.vst.plot <- ggplot(data = pca.df, aes(x = PC1, y = PC2, 
                                               color = condition,
                                               label = ID)) +
  geom_point(size = 2.5) +  
  geom_text_repel(aes(label=ID), 
                  show.legend = FALSE,
                  box.padding = unit(0.25, "lines"),
                  #fill = "white",
                  size = 4) +
  xlab(paste0(pca.x,": ", round(percentVar[pca.x] * 100), "% variance")) + 
  ylab(paste0(pca.y,": ", round(percentVar[pca.y] * 100), "% variance")) + 
  scale_color_manual(values = cond.colors) + 
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(), 
        #plot.title = element_text(colour = marks.cols.list),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 5 * .pt),
        legend.text = element_text(size = 6 * .pt),
        panel.spacing = unit(0.3, "lines"),
        strip.text.x = element_text(size = 5 * .pt),
        strip.text.y = element_text(size = 5 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3"),
        aspect.ratio = 0.8)

ggsave(filename = file.path(plots.dir,"FigS9B_RNAseq_PCA.pdf"), 
       plot = pca.vst.plot, device = cairo_pdf,
       width = 14, height = 12)



sampleDists <- dist(t(mcm2.vst))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
annot.colors <- list(condition = cond.colors)

dev.off()
cairo_pdf(filename = file.path(plots.dir, "FigS9A_RNAseq_heatmap_euclidean.pdf"),
          width = 18, height = 16)

phout <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, 
         fontsize_row = 12,
         fontsize_col = 12,
         annotation_row = samples.annot.df[,"condition",drop=FALSE],
         annotation_colors = annot.colors)#,
print(phout)
dev.off()


tmp.df <- filter(mcm2.res.df, !is.na(padj))
tmp.df$type <- "NS"
tmp.df[abs(tmp.df$log2FoldChange) >= log2FC.cutoff & tmp.df$padj < padj.cutoff, ]$type <- "DE"
tmp.df$type <- factor(tmp.df$type, levels = c("DE", "NS"))

xl <- max(abs(tmp.df$log2FoldChange))
rnaseq.volc.plot <- ggplot(tmp.df, aes(x = log2FoldChange, 
                                           y = -log10(padj),
                                           color = type)) +
  geom_vline(xintercept = log2FC.cutoff, color = "darkgray", linetype = "dashed", size = 1) + 
  geom_vline(xintercept = -log2FC.cutoff, color = "darkgray", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(padj.cutoff), color = "darkgray", linetype = "dashed", size = 1) +
  xlim(-xl,xl) +
  ylab(TeX("$-log_{10}(FDR)$")) + xlab(TeX("$log_{2}(FoldChange)$")) +
  geom_point() +
  scale_color_manual(values = c("DE" = "brown", "NS" = "gray")) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.position = "none",
        legend.text = element_text(size = 6 * .pt),
        panel.spacing = unit(0.3, "lines"),
        strip.text.x = element_text(size = 6 * .pt),
        strip.text.y = element_text(size = 6 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3"),
        aspect.ratio = 0.8)
dev.new()
ggsave(filename = file.path(plots.dir, "Fig3A_RNAseq_volcano_genes.pdf"), 
       plot = rnaseq.volc.plot, device = cairo_pdf,
       width = 14, height = 12)


mcm2.counts.df <- assay(mcm2.dds)[to.keep,]

mcm2.counts.corrected.df <- ComBat_seq(mcm2.counts.df, batch = mcm2.dds$run,
                                       group = mcm2.dds$clone)

mcm2.cpm.corrected.df <- cpm(mcm2.counts.corrected.df, 
                             lib.size = colSums(assay(mcm2.dds)),
                             log = TRUE)
wt.samps <- which(mcm2.dds$condition == "WT")

mcm2.cpm.corrected.melt.df <- melt(mcm2.cpm.corrected.df)
colnames(mcm2.cpm.corrected.melt.df) <- c("ID","sample","CPM")
mcm2.cpm.corrected.melt.df$clone <- 
  samples.annot.df[as.character(mcm2.cpm.corrected.melt.df$sample),]$clone

mcm2.cpm.corrected.mean.df <- mcm2.cpm.corrected.melt.df %>% group_by(ID, clone) %>%
  summarise(CPM = mean(CPM, na.rm=T)) %>% 
  mutate(sample = paste0("RNA_",clone)) %>%
  pivot_wider(id_cols = ID, names_from = sample, values_from = CPM) %>%
  column_to_rownames(var = "ID") %>%
  as.data.frame()

wt.clones <- paste0("RNA_",unique(as.character(samples.annot.df[samples.annot.df$condition == "WT",]$clone)))
wt.cln.smps <- match(wt.clones, colnames(mcm2.cpm.corrected.mean.df))

mcm2.zscores.corrected.df <- t(apply(mcm2.cpm.corrected.df, 1, FUN = function(x) {
  wt.mean <- mean(x[wt.samps])
  wt.sd <- sd(x[wt.samps])
  return((x - wt.mean) / wt.sd)
}))


sel.samps <- filter(samples.annot.df, replicate %in% c("r1","r2"))$ID

tmp.annot.df <- filter(samples.annot.df, replicate %in% c("r1","r2"))
tmp.cpm.corrected.df <- mcm2.cpm.corrected.df[,tmp.annot.df$ID]

wt2.samps <- which(tmp.annot.df$condition == "WT")

mcm2.zscores2.corrected.df <- t(apply(tmp.cpm.corrected.df, 1, FUN = function(x) {
  wt.mean <- mean(x[wt2.samps])
  wt.sd <- sd(x[wt2.samps])
  return((x - wt.mean) / wt.sd)
}))

  
mcm2.zscores.corrected.mean.df <- t(apply(mcm2.cpm.corrected.mean.df, 1, FUN = function(x) {
  wt.mean <- mean(x[wt.cln.smps])
  wt.sd <- sd(x[wt.cln.smps])
  return((x - wt.mean) / wt.sd)
}))

mcm2.zscores.corrected.mean.df <- t(scale(t(mcm2.cpm.corrected.mean.df)))

color.list <- list(condition = cond.colors)

fc.lim <- 3 
breakList <- seq(-fc.lim, fc.lim, by = 0.1)
pal <- maPalette(low="blue", high="red",mid="white")

hm.genes.df <- read.csv(file.path(rnaseq.data.dir,"geneList_Heatmap.csv"),
                          header = TRUE, sep = "\t")
rownames(hm.genes.df) <- hm.genes.df$gene


hm.annot.df <- filter(gene.annot.df, external_gene_name %in%  hm.genes.df$gene)
rownames(hm.annot.df) <- hm.annot.df$external_gene_name
hm.annot.df$group <- factor(hm.genes.df[hm.annot.df$external_gene_name,]$type)

clone.labels2 <- clone.labels
names(clone.labels2) <- paste0("RNA_", names(clone.labels))

tmp.mean.annot.df <- samples.annot.df %>% group_by(clone, condition) %>%
  tally() %>% mutate(ID = paste0("RNA_", clone)) %>% 
  mutate(label = clone.labels[as.character(clone)]) %>%
  column_to_rownames(var = "label") %>%
  as.data.frame()


tmp.zscores2.df <- mcm2.zscores.corrected.mean.df[hm.annot.df$ensembl_gene_id,]

gene.grp.col <- brewer.pal(length(levels(hm.annot.df$group)),
                           "Set2")
names(gene.grp.col) <- levels(hm.annot.df$group)

color.list <- list(condition = cond.colors,
                   group = gene.grp.col)

rownames(tmp.zscores2.df) <- hm.annot.df$external_gene_name
colnames(tmp.zscores2.df) <- clone.labels2[colnames(tmp.zscores2.df)]


hm.annot.df <- hm.annot.df %>% arrange(group) %>% as.data.frame()

cairo_pdf(filename = file.path(plots.dir, "FigS9D_RNAseq_heatmap_selected_genes.pdf"),
          width = 18, height = 18)
phout <- pheatmap(tmp.zscores2.df[rownames(hm.annot.df),],
         annotation_col = tmp.mean.annot.df[,c("condition"), drop = F],
         annotation_row = hm.annot.df[,c("group"), drop = F],
         cluster_cols = TRUE,
         cluster_rows = FALSE,
         clutering_method = "ward.D2",
         annotation_colors = color.list,
         color = colorRampPalette(pal)(length(breakList)),
         breaks = breakList,
         fontsize_row = 12,
         fontsize_col = 12,
         cellwidth = 24,
         cellheight = 16)

print(phout)
dev.off()



design <- model.matrix(~ run + condition, samples.annot.df)

genes.gr <- genes(gencode.txdb)
genes.gr <- subset(genes.gr, seqnames %in% chroms)

to.keep <- filterByExpr(as.matrix(gene.counts.df), design, min.count = 1)
expr.genes <- names(to.keep)[to.keep]

gene.counts.filt.df <- gene.counts.df[to.keep,]

nrst.gene <- nearest(genes.gr, RT.gr)
names(genes.gr) <- substr(names(genes.gr),1,18)

up.genes <- filter(mcm2.res.df, log2FoldChange >= log2FC.cutoff, 
                   padj < padj.cutoff)$ensembl_gene_id
down.genes <- filter(mcm2.res.df, log2FoldChange <= -log2FC.cutoff, 
                     padj < padj.cutoff)$ensembl_gene_id
genes.gr$RT <- RT.gr$RT.class[nrst.gene]
genes.gr$DE <- "No"
genes.gr[up.genes]$DE <- "Up"
genes.gr[down.genes]$DE <- "Down"

genes.gr <- genes.gr[names(genes.gr) %in% expr.genes]

gene.rt.df <- c()
for (dt in c("Up","Down")) {
  dt.genes <- names(subset(genes.gr, DE == dt))
  for (rt in c("early","mid-early","mid-late","late")) {
    rt.genes <- names(subset(genes.gr, RT == rt))
    dt.and.rt <- intersect(dt.genes, rt.genes)
    dt.not.rt <- setdiff(dt.genes, rt.genes)
    rt.not.dt <- setdiff(rt.genes, dt.genes)
    the.rest <- setdiff(names(genes.gr), union(dt.genes, rt.genes))
    tb <- matrix(c(length(dt.and.rt), length(rt.not.dt), 
                   length(dt.not.rt), length(the.rest)), nrow = 2)
    ft <- fisher.test(tb)
    tmp.df <- data.frame(DE = dt, RT = rt, pvalue = ft$p.value, odds = ft$estimate,
                         lower.conf = ft$conf.int[1], upper.conf = ft$conf.int[2])
    tmp.df$n <- length(dt.and.rt)
    gene.rt.df <- rbind(gene.rt.df, tmp.df)
  }
}
gene.rt.df$DE <- factor(gene.rt.df$DE, levels = c("Down","Up"))
gene.rt.df$RT <- factor(gene.rt.df$RT, levels = c("early","mid-early","mid-late","late"))

gene.rt.df$label <- trimws(paste(gene.rt.df$n, pval2signf(gene.rt.df$pvalue)))
gene.rt.df$log2odds <- log2(gene.rt.df$odds)
gene.rt.df$padj <- p.adjust(gene.rt.df$pvalue)
gene.rt.df$log2odds[gene.rt.df$padj > 0.001] <- NA

fig.rt.genes.plt <- ggplot(gene.rt.df,aes(x = RT, y = DE, fill = log2odds, 
                      label = n)) +
  geom_tile() + 
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-2,2)) +
  ylab("DE genes") + xlab("Replication Timing") +
  labs(fill = "Log-odds") +
  geom_text(aes(label=n), size=12) +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  theme(panel.spacing = unit(0.3, "lines"),
        legend.title = element_text(size = 10 * .pt),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12 * .pt),
        axis.text = element_text(size = 10 * .pt),
        legend.text = element_text(size = 7 * .pt))

ggsave(filename = file.path(plots.dir, "FigS10C_RNAseq_de_genes_RT.pdf"), 
       plot = fig.rt.genes.plt, device = cairo_pdf,
       width = 15, height = 8)

### ================= MCM2-2A Repeat analyses ==================================


te.counts.df <- read.csv(file.path(rnaseq.data.dir,"MCM2_TEsubf_count_matrix.tsv.gz"),
                           header = TRUE, sep = "\t")

rownames(te.counts.df) <- te.counts.df$gene.id
te.counts.df <- te.counts.df[,-1]


mcm2.te.dds <- DESeqDataSetFromMatrix(countData = as.matrix(te.counts.df),
                                   colData = samples.annot.df[colnames(te.counts.df),],
                                   rowData = gene.annot.df[rownames(te.counts.df),],
                                   design = ~ run + condition)

mcm2.te.ddq <- DESeq(mcm2.te.dds)
mcm2.te.res <- results(mcm2.te.ddq, name = "condition_MCM2_2A_vs_WT")
mcm2.te.res.df <- as.data.frame(mcm2.te.res)
mcm2.te.res.df <- cbind(mcm2.te.res.df, as.data.frame(rowData(mcm2.te.dds)))


tmp.df <- filter(mcm2.te.res.df, !is.na(padj), !is.na(Repeat.Subfamily))
tmp.df$type <- "No.diff"
tmp.df[tmp.df$log2FoldChange >= log2FC.cutoff & tmp.df$padj < padj.cutoff, ]$type <- "DE.up"
tmp.df$type <- factor(tmp.df$type)
tmp.df$rep.label <- sapply(tmp.df$ID, FUN = function(x) {
  return(unlist(strsplit(x, ":"))[1])
})
tmp.df$label <- ""
tmp.df[tmp.df$type == "DE.up",]$label <- tmp.df[tmp.df$type == "DE.up",]$rep.label
tmp.df$family.color <- as.character(tmp.df$Repeat.Family)
fms <- c(unique(tmp.df$family.color),"NS")
tmp.df[tmp.df$label == "",]$family.color <- "NS"
tmp.df$family.color <- factor(tmp.df$family.color, levels = c("ERV1","ERVK","ERVL","Satellite","NS"))

rep.de.fam.cols <- brewer.pal(8, "Dark2")[c(1,7,3,4,8)]
names(rep.de.fam.cols) <- levels(tmp.df$family.color)

xl <- max(abs(tmp.df$log2FoldChange))
xl <- 3
tb.cnts <- table(tmp.df$type)
rnaseq.rep.mcm2.plot <- ggplot(tmp.df, aes(x = log2FoldChange, 
                                           y = -log10(padj),
                                           color = family.color,
                                           label = label)) +
  geom_vline(xintercept = log2FC.cutoff, color = "darkgray", linetype = "dashed", size = 1) + 
  geom_vline(xintercept = -log2FC.cutoff, color = "darkgray", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(padj.cutoff), color = "darkgray", linetype = "dashed", size = 1) +
  xlim(-xl,xl) +
  ylab(TeX("$-log_{10}(FDR)$")) + xlab(TeX("$log_{2}(FoldChange)$")) +
  geom_label_repel(size = 6, show.legend = FALSE) +
  geom_point() + 
  scale_color_manual(values = rep.de.fam.cols) +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  guides(color=guide_legend(title="Repeat family")) +
  theme(axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        panel.spacing = unit(0.3, "lines"),
        strip.text.x = element_text(size = 6 * .pt),
        strip.text.y = element_text(size = 6 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3"),
        aspect.ratio = 0.8)
ggsave(filename = file.path(plots.dir, "Fig3B_RNAseq_repeats_volcano.pdf"), 
       width = 18, height = 11, 
       plot = rnaseq.rep.mcm2.plot, device = cairo_pdf)



to.keep <- rowSums(assay(mcm2.te.dds)) > 1
names(to.keep) <- rownames(mcm2.te.dds)
  
mcm2.te.counts.df <- assay(mcm2.te.dds)[to.keep,]

mcm2.te.counts.corrected.df <- ComBat_seq(mcm2.te.counts.df, 
                                          batch = mcm2.te.dds$run,
                                          group = mcm2.te.dds$clone)

mcm2.te.cpm.corrected.df <- cpm(mcm2.te.counts.corrected.df, 
                                lib.size = colSums(assay(mcm2.te.dds)),
                                log = TRUE)

wt.samps <- which(mcm2.te.dds$condition == "WT")

mcm2.te.cpm.corrected.melt.df <- melt(mcm2.te.cpm.corrected.df)
colnames(mcm2.te.cpm.corrected.melt.df) <- c("ID","sample","CPM")
mcm2.te.cpm.corrected.melt.df$clone <- 
  samples.annot.df[as.character(mcm2.te.cpm.corrected.melt.df$sample),]$clone

mcm2.te.cpm.corrected.mean.df <- mcm2.te.cpm.corrected.melt.df %>% group_by(ID, clone) %>%
  summarise(CPM = mean(CPM, na.rm=T)) %>% 
  mutate(sample = paste0("RNA_",clone)) %>%
  pivot_wider(id_cols = ID, names_from = sample, values_from = CPM) %>%
  column_to_rownames(var = "ID") %>%
  as.data.frame()

wt.clones <- paste0("RNA_",unique(as.character(samples.annot.df[samples.annot.df$condition == "WT",]$clone)))
wt.cln.smps <- match(wt.clones, colnames(mcm2.te.cpm.corrected.mean.df))
sel.samps <- filter(samples.annot.df, replicate %in% c("r1","r2"))$ID
wt2.samps <- which(tmp.annot.df$condition == "WT")

mcm2.te.zscores.corrected.mean.df <- t(apply(mcm2.te.cpm.corrected.mean.df, 1, FUN = function(x) {
  wt.mean <- mean(x[wt.cln.smps])
  wt.sd <- sd(x[wt.cln.smps])
  return((x - wt.mean) / wt.sd)
}))

mcm2.te.zscores.corrected.mean.df <- t(scale(t(mcm2.te.cpm.corrected.mean.df)))

de.te <- filter(mcm2.te.res.df, abs(log2FoldChange) >= log2FC.cutoff, 
                padj < padj.cutoff,
                !is.na(Repeat.Subfamily))$ID




tmp.te.annot.df <- as.data.frame(rowData(mcm2.te.dds))[de.te,]
tmp.zscores.df <- mcm2.te.zscores.corrected.mean.df[tmp.te.annot.df$ID,]

rownames(tmp.te.annot.df) <- tmp.te.annot.df$Repeat.Subfamily
rownames(tmp.zscores.df) <- tmp.te.annot.df$Repeat.Subfamily
colnames(tmp.zscores.df) <- clone.labels2[colnames(tmp.zscores.df)]


tmp.te.annot.df$Family <- tmp.te.annot.df$Repeat.Family
tmp.te.annot.df$Family <- factor(tmp.te.annot.df$Family, levels = c("ERV1","ERVK",
                                                                    "ERVL","Satellite"))

tmp.te.annot.df <- tmp.te.annot.df %>% arrange(Family) %>% as.data.frame()

gene.grp.col <- brewer.pal(length(levels(alva.annot.df$group)),
                           "Set2")
names(gene.grp.col) <- levels(alva.annot.df$group)

te.color.list <- list(condition = cond.colors,
                      Family = rep.de.fam.cols[levels(tmp.te.annot.df$Family)])


tmp.zscores.df <- tmp.zscores.df[rownames(tmp.te.annot.df),]

cairo_pdf(filename = file.path(plots.dir, "FigS9E_RNAseq_heatmap_DE_Repeats.pdf"),
          width = 18, height = 16)

phout <- pheatmap(tmp.zscores.df,
                  annotation_col = tmp.mean.annot.df[,c("condition"), drop = F],
                  annotation_row = tmp.te.annot.df[,c("Family"), drop = F],
                  cluster_rows = FALSE,
                  clustering_method = "ward.D2",
                  clustering_distance_cols = "euclidean",
                  annotation_colors = te.color.list,
                  color = colorRampPalette(pal)(length(breakList)),
                  breaks = breakList,
                  fontsize_row = 12,
                  fontsize_col = 12,
                  cellwidth = 22,
                  cellheight = 22)

print(phout)
dev.off()


### ================= POLE4-KO datasets =========================================


pole4.gene.counts.df <- read.csv(file = file.path(rnaseq.data.dir,"POLE4_gene_count_matrix.tsv.gz"),
                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(pole4.gene.counts.df) <- pole4.gene.counts.df$gene.id
pole4.gene.counts.df <- pole4.gene.counts.df[,-1]

pole4.te.counts.df <- read.csv("data/RNA_seq/POLE4_TEsubf_count_matrix.tsv.gz",
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(pole4.te.counts.df) <- pole4.te.counts.df$gene.id
pole4.te.counts.df <- pole4.te.counts.df[,-1]

pole4.samples.df <- read.csv("data/RNA_seq/POLE4_samples.tsv",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(pole4.samples.df) <- pole4.samples.df$ID

pole4.samples.df$condition <- factor(pole4.samples.df$condition, 
                                     levels = c("WT","MCM2_2A","POLE4_KO"))

pole4.samples.df$clone <- factor(pole4.samples.df$clone, 
                                     levels = c("410","439","550","551","552"))


pole4.samples.df$replicate <- factor(pole4.samples.df$replicate)

pole4.dds <- DESeqDataSetFromMatrix(countData = pole4.gene.counts.df,
                                    rowData = gene.annot.df[rownames(pole4.gene.counts.df),],
                                    colData = pole4.samples.df[colnames(pole4.gene.counts.df),],
                                    design = ~ clone)


pole4.tmp.samples.df <- filter(pole4.samples.df,condition != "MCM2_2A")
pole4.tmp.counts.df <- pole4.te.counts.df[,pole4.tmp.samples.df$ID]
pole4.tmp.samples.df$clone <- droplevels(pole4.tmp.samples.df$clone)
pole4.tmp.samples.df$condition <- droplevels(pole4.tmp.samples.df$condition)

pole4.dds <- DESeqDataSetFromMatrix(countData = pole4.tmp.counts.df,
                                    rowData = gene.annot.df[rownames(pole4.tmp.counts.df),],
                                    colData = pole4.tmp.samples.df[colnames(pole4.tmp.counts.df),],
                                    design = ~ clone)

pole4.ddq <- DESeq(pole4.dds)
pole4.de.gene.list <- list()
pole4.results.list <- list()
for (cid in c("550","551","552")) {
  pole4.res <- results(pole4.ddq, name = paste0("clone_",cid,"_vs_410"))
  pole4.res.df <- as.data.frame(pole4.res)
  pole4.res.df <- cbind(pole4.res.df, as.data.frame(rowData(pole4.dds)))
  pole4.de.gene.list[[cid]] <- filter(pole4.res.df, abs(log2FoldChange) >= log2FC.cutoff,
                                padj < padj.cutoff)$ensembl_gene_id
  pole4.results.list[[cid]] <- pole4.res.df
}


design(pole4.dds) <- ~ condition
pole4.ddq <- DESeq(pole4.dds)
pole4.res <- results(pole4.ddq, name = "condition_POLE4_KO_vs_WT")
pole4.res.df <- as.data.frame(pole4.res)
pole4.res.df <- cbind(pole4.res.df, as.data.frame(rowData(pole4.dds)))
pole4.de.gene.list[["full"]] <- filter(pole4.res.df, abs(log2FoldChange) >= log2FC.cutoff,
                                    padj < padj.cutoff, is.na(Repeat.Class))$ID

pole4.results.list[["full"]] <- pole4.res.df


mcm2.te.res.df <- mcm2.te.res.df[rownames(pole4.res.df),]

mcm2.pole4.cor.df <- data.frame(MCM2_2A.log2FC = mcm2.te.res.df$log2FoldChange,
                                POLE4_KO.log2FC = pole4.res.df$log2FoldChange,
                                MCM2_2A.padj = mcm2.te.res.df$padj,
                                POLE4_KO.padj = pole4.res.df$padj, stringsAsFactors = FALSE)
tmp.df <- filter(mcm2.pole4.cor.df, (abs(MCM2_2A.log2FC) >= log2FC.cutoff), 
                (abs(POLE4_KO.log2FC) >= log2FC.cutoff), (!is.na(MCM2_2A.padj)) |
                  (!is.na(POLE4_KO.padj))) %>% as.data.frame()

max.x <- max(abs(tmp.df$MCM2_2A.log2FC), na.rm = TRUE)
max.y <- max(abs(tmp.df$POLE4_KO.log2FC), na.rm = TRUE)



p4mcm.cor <- cor.test(tmp.df$MCM2_2A.log2FC, tmp.df$POLE4_KO.log2FC, 
                      method = "spearman")

fig.mcm2.vs.pole4.plt <- 
  ggplot(tmp.df, aes(x = MCM2_2A.log2FC, y = POLE4_KO.log2FC)) +
  geom_smooth(method='lm', formula = y~x, color = "blue") +
  geom_point(size = 2, shape = 1) +
  geom_density2d() +
  xlim(-max.x, max.x) +
  ylim(-max.y, max.y) +
  xlab("MCM2-2A Log2(Fold Change") +
  ylab("POLE4-KO Log2(Fold Change") +
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(filename = file.path(plots.dir, "Fig3G_RNAseq_MCM2.vs.POLE4.pdf"), 
       width = 14, height = 12, 
       plot = fig.mcm2.vs.pole4.plt, device = cairo_pdf)


### ================= External datasets =========================================


public.counts.df <- read.csv(file.path(rnaseq.data.dir,"RNAseq_public_counts.tsv.gz"),
                             sep = "\t", header = TRUE)

rownames(public.counts.df) <- public.counts.df$gene.id
public.counts.df <- public.counts.df[,-1]

public.samples.df <- read.csv(file.path(rnaseq.data.dir,"RNAseq_public_samples.tsv"),
                               sep = "\t", header = TRUE)

rownames(public.samples.df) <- public.samples.df$ID

public.samples.df$condition <- gsub("-","_", public.samples.df$condition)


suv39.samples.df <- filter(public.samples.df, accession == "GSE57092")
suv39.samples.df$condition <- factor(suv39.samples.df$condition, 
                                     levels = c("WT","SUV39H_dKO"))
suv39.samples.df$replicate <- factor(suv39.samples.df$replicate)

suv39.dds <- DESeqDataSetFromMatrix(countData = public.counts.df[,suv39.samples.df$ID],
                                    rowData = gene.annot.df[rownames(public.counts.df),],
                                    colData = suv39.samples.df,
                                    design = ~ condition)
suv39.ddq <- DESeq(suv39.dds)
suv39.res.df <- as.data.frame(results(suv39.ddq))
suv39.res.df <- cbind(suv39.res.df, gene.annot.df[rownames(suv39.res.df),])
suv39.res.df$study <- "SUV39H_dKO"


suz12.samples.df <- filter(public.samples.df, accession == "GSE127804")
suz12.samples.df$condition <- factor(suz12.samples.df$condition, 
                                     levels = c("WT","SUZ12_KO"))
suz12.samples.df$replicate <- factor(suz12.samples.df$replicate)

suz12.dds <- DESeqDataSetFromMatrix(countData = public.counts.df[,suz12.samples.df$ID],
                                    rowData = gene.annot.df[rownames(public.counts.df),],
                                    colData = suz12.samples.df,
                                    design = ~ condition)
suz12.ddq <- DESeq(suz12.dds)
suz12.res.df <- as.data.frame(results(suz12.ddq))
suz12.res.df <- cbind(suz12.res.df, gene.annot.df[rownames(suz12.res.df),])
suz12.res.df$study <- "SUZ12_KO"


setdb1.samples.df <- filter(public.samples.df, accession == "PRJNA544540")
setdb1.samples.df$condition <- factor(setdb1.samples.df$condition, 
                                     levels = c("SETDB1_KO_no4OHT","SETDB1_KO_4OHT"))
setdb1.samples.df$replicate <- factor(setdb1.samples.df$replicate)

setdb1.dds <- DESeqDataSetFromMatrix(countData = public.counts.df[,setdb1.samples.df$ID],
                                    rowData = gene.annot.df[rownames(public.counts.df),],
                                    colData = setdb1.samples.df,
                                    design = ~ condition)
setdb1.ddq <- DESeq(setdb1.dds)
setdb1.res.df <- as.data.frame(results(setdb1.ddq))
setdb1.res.df <- cbind(setdb1.res.df, gene.annot.df[rownames(setdb1.res.df),])
setdb1.res.df$study <- "SETDB1_KO"



gene.de.up.list <- list("MCM2-2A" = mcm2.gene.de.list$Up,
                        "SUV39H-dKO" = filter(suv39.res.df, log2FoldChange >= log2FC.cutoff,
                                              padj < padj.cutoff, is.na(Repeat.Class))$ID,
                        "SETDB1-KO" = filter(setdb1.res.df, log2FoldChange >= log2FC.cutoff,
                                             padj < padj.cutoff, is.na(Repeat.Class))$ID,
                        "SUZ12-KO" =  filter(suz12.res.df, log2FoldChange >= log2FC.cutoff,
                                             padj < padj.cutoff, is.na(Repeat.Class))$ID)

te.de.up.list <- list("MCM2-2A" = filter(mcm2.te.res.df, log2FoldChange >= log2FC.cutoff,
                                           padj < padj.cutoff, !is.na(Repeat.Class))$ID,
                        "SUV39H-dKO" = filter(suv39.res.df, log2FoldChange >= log2FC.cutoff,
                                              padj < padj.cutoff, !is.na(Repeat.Class))$ID,
                        "SETDB1-KO" = filter(setdb1.res.df, log2FoldChange >= log2FC.cutoff,
                                             padj < padj.cutoff, !is.na(Repeat.Class))$ID,
                        "SUZ12-KO" =  filter(suz12.res.df, log2FoldChange >= log2FC.cutoff,
                                             padj < padj.cutoff, !is.na(Repeat.Class))$ID)


gene.k9f.plt <- ggvenn(gene.de.up.list[c("MCM2-2A","SUV39H-dKO","SETDB1-KO")],
       text_size = 6) +
  theme_void(base_size = 18, base_family = "Helvetica")

ggsave(filename = file.path(plots.dir, "FigS11D_RNAseq_gene_up_k9fact_venn.pdf"), 
       width = 14, height = 12, 
       plot = gene.k9f.plt, device = cairo_pdf)


gene.k27f.plt <- ggvenn(gene.de.up.list[c("MCM2-2A","SUZ12-KO")],
                       text_size = 6) +
  theme_void(base_size = 18, base_family = "Helvetica")

ggsave(filename = file.path(plots.dir, "FigS11E_RNAseq_gene_up_k27fact_venn.pdf"), 
       width = 14, height = 12, 
       plot = gene.k27f.plt, device = cairo_pdf)

te.k9f.plt <- ggvenn(te.de.up.list[c("MCM2-2A","SUV39H-dKO","SETDB1-KO")],
                       text_size = 6) +
  theme_void(base_size = 18, base_family = "Helvetica")

ggsave(filename = file.path(plots.dir, "FigS11F_RNAseq_TE_up_k9fact_venn.pdf"), 
       width = 14, height = 12, 
       plot = te.k9f.plt, device = cairo_pdf)


