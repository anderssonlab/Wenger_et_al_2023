library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(stringi)
library(ChIPseeker)
library(foreach)
library(doParallel)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(scales)
library(reshape2)
library(bamsignals)
library(csaw)
library(BiocParallel)
library(clusterProfiler)
library(limma)
library(tidyverse)
library(Rsamtools)
library(ChIPpeakAnno)
library(RColorBrewer)
library(edgeR)

library(ggplot2)
library(ggforce)
library(ggpmisc)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(UpSetR)
library(ecoflux)
library(seqsetvis)
library(pheatmap)
library(iterators)
library(parallel)
library(foreach)
library(ggh4x)
library(coin)
library(rstatix)
library(openxlsx)

## ======================= Load annotations ====================================
wd.dir <- "Wenger_et_al_2022"
setwd(wd.dir)

data.dir <- "data"
plots.dir <- "plots"

chip.data.dir <- file.path(data.dir,"ChIP_seq")
external.data.dir <- "data/external"

## Load process external data and functions
load(file.path(external.data.dir, "Processed_external.RData"))
source(file.path(wd.dir,"utils/utils.R"))
num.cores <- detectCores() / 2


annot.version <- "gencode.vM23"
assembly <- "mm10"

chrom.sizes.file <- file.path(data.dir, paste0(assembly, ".chrom.sizes"))
chrom.sizes.df <- read.csv2(file = chrom.sizes.file, header = FALSE, col.names = c("chrom", "length"), sep = "\t")

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

## =================== Define clone labels =====================================
reps <- c("r1" = "r1", "r2" = "r2", "r3" = "r1", "r4" = "r1", "r5" = "r2")

clone.labs <- c("WT","MCM2-2A #1","MCM2-2A #2", "MCM2-2A #3", "MCM2-2A #4",
                "MCM2-Y90A",
                "MCM2-R #1", "MCM2-R #2","POLE4-KO #1", "POLE4-KO #2")

names(clone.labs) <- c("410","439","442","438", 
                       "421","423",
                       "586","588","551","552")  

clone.conds <- c("WT", "MCM2-2A", "MCM2-2A","MCM2-2A", 
                 "MCM2-2A","MCM2-Y90A","MCM2-R","MCM2-R",
                 "POLE4-KO", "POLE4-KO")

names(clone.conds) <- names(clone.labs)

replicate.labs <- c("replicate #1", "replicate #2", "replicate #1","replicate #1", 
                    "replicate #2")
names(replicate.labs) <- c("r1", "r2", "r3","r4","r5")


tmp.df <- crossing(names(clone.labs), names(replicate.labs)) %>% as.data.frame()
colnames(tmp.df) <- c("clone","replicate")
tmp.df$sample <- paste0(tmp.df[,1], "_", tmp.df[,2])
tmp.df$label <- paste0(clone.labs[tmp.df[,1]], "\n", replicate.labs[tmp.df[,2]])

tmp.df$clone <- factor(tmp.df$clone, levels =  names(clone.labs))
tmp.df$replicate <- factor(tmp.df$replicate, 
                           levels =  c("r1","r3","r4","r5","r2"))

tmp.df <- tmp.df %>% arrange(clone, replicate) %>% as.data.frame()
samp.order <- tmp.df$sample

samp.label <- tmp.df$label
names(samp.label) <- tmp.df$sample

annot.col <- brewer.pal(7, "Set2")
names(annot.col) <- c("Promoter","exon","Distal Intergenic",
                      "intron","5' UTR","3' UTR","Downstream <3kb")


## ===================  Load external data =====================================
chroms <- c(paste("chr",1:19,sep=""), "chrX", "chrY")
spike.chroms <- c("chr2L_spike", "chr2LHet_spike","chr2R_spike", "chr2RHet_spike", 
                  "chr3L_spike", "chr3LHet_spike", "chr3R_spike", "chr3RHet_spike", 
                  "chr4_spike", "chrM_spike", "chrU_spike", "chrUextra_spike",
                  "chrX_spike", "chrXHet_spike","chrYHet_spike")

chip.external.data.dir <- "data/ChIP_seq"
plots.dir <- "plots"

### Load processed peak ranges
load(file.path(chip.data.dir,"ChIP_processed_peaks.RData"))

mark.list <- c("H3K27ac","H3K4me3","H3K27me3","H3K9me3")
clone.list <- c("410","421","439","442")

pkset.names <- c("410" = "Lost_3", "410///421" = "Lost_2",
                 "410///439" = "Lost_2", "410///442" = "Lost_2",
                 "410///421///439" = "Lost_1", 
                 "410///421///442" = "Lost_1",
                 "410///439///442" = "Lost_1",
                 "410///421///439///442" = "Kept",
                 "421" = "Gained_1", "439" = "Gained_1","442" = "Gained_1",
                 "421///439" = "Gained_2", "421///442" = "Gained_2",
                 "439///442" = "Gained_2", "421///439///442" = "Gained_3")


external.marks <- c("H2AK119ub1","H3K36me3","H3K27me1",
                    "H3K27me2", "H3K9me2", "H3K36me2","H3K4me1")
our.marks <- c("H3K4me3","H3K27ac","H3K9me3","H3K27me3")
peak.set.lvls <- c("Lost_3","Lost_2","Lost_1","Kept","Gained_1","Gained_2","Gained_3")
peak.set.lvls3 <- c("Lost","Kept","Gained")

### Peak count analyses
all.pk.mark.list <- list()

all.peak.sets.df <- c()
for (mk in mark.list) {
  print(mk)
  mids <- paste0(mk, "_",clone.list)
  pk.ovrl <- findOverlapsOfPeaks(chip.peaks.list[mids], minoverlap = 1,
                                 connectedPeaks = "merge")
  
  names(pk.ovrl$peaklist) <- gsub(paste0(mk,"_"),"",names(pk.ovrl$peaklist))
  mk.gr <- GRanges()
  for (pkn in names(pk.ovrl$peaklist)) {
    tmp.gr <- pk.ovrl$peaklist[[pkn]]
    tmp.gr$peak.7class <- pkset.names[pkn]
    mk.gr <- c(mk.gr, tmp.gr)
  }
  mk.gr <- annotateRegions(mk.gr, gencode.txdb, 
                           tss.region = c(-3000,3000),
                           Enh.gr, SEnh.gr, cpg.islands.gr,
                           all_CpG_sites, cpg.meth.gr)
  mk.gr$mCpG.island <- mk.gr$mCpG.fraction & mk.gr$CpG.island
  mk.gr$uCpG.island <- !mk.gr$mCpG.fraction & mk.gr$CpG.island
  
  for (ext.mark in external.marks) {
    elementMetadata(mk.gr)[[ext.mark]] <- overlapsAny(mk.gr, chip.external.peaks.list[[ext.mark]],
                                                      minoverlap = 1)
  }
  
  for (other.mark in setdiff(our.marks,mk)) {
    elementMetadata(mk.gr)[[other.mark]] <- overlapsAny(mk.gr, 
                                                        chip.peaks.list[[paste0(other.mark,"_410")]],
                                                        minoverlap = 1)
  }
  
  #all.mks <- c(external.marks, setdiff(our.marks, mk), "CpG.island","mCpG.island")
  all.mks <- c(external.marks, setdiff(our.marks, mk), "mCpG.island")
  elementMetadata(mk.gr)[["None"]] <- rowSums(as.data.frame(mk.gr)[,all.mks]) == 0
  mk.gr$peak.7class <- factor(mk.gr$peak.7class, levels = peak.set.lvls)
  
  rt.dist <- as.data.frame(distanceToNearest(mk.gr, RT.gr))
  
  elementMetadata(mk.gr)[["RT.class"]] <- RT.gr[rt.dist$subjectHits]$RT.class
  peak.sets.df <- as.data.frame(mk.gr)[,c("annotation_cur","enhancer.class",
                                          "peak.7class","RT.class")]
  peak.sets.df$mark <- mk
  peak.sets.df$width <- width(mk.gr)
  peak.sets.df$GC.content <- mk.gr$GC.cont
  peak.sets.df$distToIZnorm <- mk.gr$distToIZnorm
  peak.sets.df$distToTSS <- mk.gr$distToTSS
  
  all.peak.sets.df <- rbind(all.peak.sets.df, peak.sets.df)
  all.pk.mark.list[[mk]] <- mk.gr
}



peak.set.cols <- c(brewer.pal(9,"Blues")[c(7,5,3)],"gray",
                   brewer.pal(9,"Oranges")[c(3,5,7)])

names(peak.set.cols) <- peak.set.lvls


peak.set.cols3 <- peak.set.cols[c(1,4,7)]
names(peak.set.cols3) <- peak.set.lvls3

peak.set.cols.list <- list("7class" = peak.set.cols,
                           "3class" = peak.set.cols3)


all.peak.sets.df$peak.group <- all.peak.sets.df$peak.7class

peak.sets.mean.df <- all.peak.sets.df %>% group_by(mark,peak.group) %>%
  tally() %>% mutate(fraction = n / sum(n), percentage = (n / sum(n)) * 100) %>%
  as.data.frame()

peak.sets.mean.df$label <- paste0(sprintf("%.1f", peak.sets.mean.df$percentage),"%")

all.peak.sets.mean.df <- c()
for (mk in unique(peak.sets.mean.df$mark)) {
  df2 <- filter(peak.sets.mean.df, mark == mk) %>% 
    mutate(csum = rev(cumsum(rev(percentage))), 
           pos = percentage/2 + lead(csum, 1)) %>%
    mutate(pos = if_else(is.na(pos), percentage/2, pos)) %>% 
    as.data.frame()
  all.peak.sets.mean.df <- rbind(all.peak.sets.mean.df, df2)
}


chip.fig2A.plt <- ggplot(all.peak.sets.mean.df, 
                      aes(x = peak.group , y = percentage, fill = peak.group, label = n)) +
  geom_bar(width = 0.5, position = "dodge", stat = "identity") +
  scale_fill_manual(values = peak.set.cols.list$`7class`) +  
  geom_text(aes(label = n), 
            vjust = -0.15,
            size = 5) +
  geom_point(aes(y = percentage + 4.5), color = "white", alpha = 0) +
  ylab("Percentage (%) of peaks") +
  xlab("") +
  facet_wrap(~mark, nrow = 2) +#, scales = "free_y") +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 5 * .pt),
        legend.text = element_text(size = 6 * .pt),
        strip.text.x = element_text(size = 6 * .pt),
        strip.text.y = element_text(size = 6 * .pt),
        strip.background = element_rect(colour="black", fill="lightgray"),
        panel.spacing = unit(0.3, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5 * .pt),
        aspect.ratio = 0.6)

ggsave(filename = file.path(plots.dir,"Fig2A_ChIP_peak_counts.pdf"),
       plot = chip.fig2A.plt, 
       width = 14, height = 10,
       device = cairo_pdf)



### Do upset plots


cids <- c("410","421","439","442")
clone.mark.ovl.list <- mclapply(mark.list, mc.cores = min(4,num.cores), 
                                FUN = function(mk) {
  pks <- paste(mk, cids, sep = "_")
  tmp.list <- chip.peaks.list[pks]
  pk.ovrl <- findOverlapsOfPeaks(tmp.list,connectedPeaks = "merge", 
                                 minoverlap = 1)
})
names(clone.mark.ovl.list) <- mark.list

for (mk in mark.list) {
  print(mk)
  tmp.mat <- matrix(clone.mark.ovl.list[[mk]]$venn_cnt, 
                    nrow = nrow(clone.mark.ovl.list[[mk]]$venn_cnt))
  
  colnames(tmp.mat) <- colnames(clone.mark.ovl.list[[mk]]$venn_cnt)
  
  upset.mat <- tmp.mat[unlist(sapply(1:nrow(tmp.mat),FUN = function(x){
    if (tmp.mat[x,5] > 0) {
      return(rep(x,tmp.mat[x,5]))
    }
  })),1:4]
  
  colnames(upset.mat) <- cids
  rownames(upset.mat) <- paste0("Peak_",1:nrow(upset.mat))
  
  ovl.pk.list <- list()
  for (cid in colnames(upset.mat)) {
    ovl.pk.list[[cid]] <- rownames(upset.mat)[upset.mat[,cid] == 1]
  }
  names(ovl.pk.list) <- clone.labs[names(ovl.pk.list)]
  upset.mat <- fromList(ovl.pk.list)
  inter.list <- list(list(clone.labs["410"]),
                     list(clone.labs["410"],clone.labs["421"]),
                     list(clone.labs["410"],clone.labs["439"]),
                     list(clone.labs["410"],clone.labs["442"]),
                     list(clone.labs["410"],clone.labs["421"],clone.labs["439"]),
                     list(clone.labs["410"],clone.labs["421"],clone.labs["442"]),
                     list(clone.labs["410"],clone.labs["439"],clone.labs["442"]),
                     list(clone.labs["410"],clone.labs["421"],clone.labs["439"],clone.labs["442"]),
                     list(clone.labs["421"]),list(clone.labs["439"]),list(clone.labs["442"]),
                     list(clone.labs["421"],clone.labs["439"]),list(clone.labs["439"],clone.labs["442"]),list(clone.labs["421"],clone.labs["442"]),
                     list(clone.labs["421"],clone.labs["439"],clone.labs["442"]))
  
  cairo_pdf(filename = file.path(plots.dir,
                                 paste0("FigS5A_ChIP_",mk,"_upset.pdf")),
            width = 14, height = 10)
  
  up.plot <- upset(upset.mat, sets = colnames(upset.mat),
                   point.size = 2.5, line.size = 0.9,
                   keep.order = TRUE, 
                   order.by = "freq",
                   text.scale = 2,
                   mainbar.y.label = "Number of Overlaps",
                   sets.x.label = paste0(mk, " Peaks"),
                   sets.bar.color = c("blue","orange","orange","orange"),
                   empty.intersections = "on",
                   #intersections = inter.list,
                   queries = list(list(query = intersects, 
                                       params = list(clone.labs["410"]), color = peak.set.cols["Lost_3"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["410"],clone.labs["421"]), color = peak.set.cols["Lost_2"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["410"],clone.labs["439"]), color = peak.set.cols["Lost_2"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["410"],clone.labs["442"]), color = peak.set.cols["Lost_2"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["410"],clone.labs["421"],clone.labs["439"]), color = peak.set.cols["Lost_1"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["410"],clone.labs["421"],clone.labs["442"]), color = peak.set.cols["Lost_1"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["410"],clone.labs["439"],clone.labs["442"]), color = peak.set.cols["Lost_1"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["410"],clone.labs["421"],clone.labs["439"],clone.labs["442"]), color = peak.set.cols["Kept"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["421"]), color = peak.set.cols["Gained_1"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["439"]), color = peak.set.cols["Gained_1"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["442"]), color = peak.set.cols["Gained_1"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["421"],clone.labs["439"]), color = peak.set.cols["Gained_2"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["421"],clone.labs["442"]), color = peak.set.cols["Gained_2"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["439"],clone.labs["442"]), color = peak.set.cols["Gained_2"], active = T),
                                  list(query = intersects, 
                                       params = list(clone.labs["421"],clone.labs["439"],clone.labs["442"]), color = peak.set.cols["Gained_3"], active = T)))
  print(up.plot)
  dev.off()
}



### Chromatin state analysis of peaks



chrom.marks <- c("mCpG.island",
                 "H3K9me3","H3K9me2",
                 "H3K27me3","H3K27me2","H3K27me1",
                 "H2AK119ub1","H3K36me3","H3K36me2",
                 "H3K4me3","H3K4me1","H3K27ac","None")
ft.chr.df <- c()
ft.annot.df <- c()

for (mk in names(all.pk.mark.list)) {
  tmp.gr <- all.pk.mark.list[[mk]]
  tmp.gr$peak.group <- tmp.gr$peak.7class
  tmp2.chr.df <- c()
  tmp2.annot.df <- c()
  for (pk.grp in unique(tmp.gr$peak.group)) {
    for (annot in unique(tmp.gr$annotation_cur)) {
      inft.inst <- length(subset(tmp.gr, (annotation_cur == annot) &
                                   (peak.group == pk.grp)))
      inft.onst <- length(subset(tmp.gr, (annotation_cur == annot) &
                                   (peak.group != pk.grp)))
      onft.inst <- length(subset(tmp.gr, (annotation_cur != annot) &
                                   (peak.group == pk.grp)))
      
      onft.onst <- length(subset(tmp.gr, (annotation_cur == annot) |
                                   (peak.group == pk.grp)))
      
      am <- matrix(c(inft.inst,
                     inft.onst,
                     onft.inst,
                     length(tmp.gr) - onft.onst),
                   byrow = TRUE, ncol = 2)
      ft <- fisher.test(am)
      tmp.annot.df <- data.frame(annotation = annot,
                                 peak.group = pk.grp, 
                                 mark = mk, 
                                 pvalue = ft$p.value,
                                 odds.ratio = ft$estimate,
                                 log2odds = log2(ft$estimate),
                                 conf.95.low = ft$conf.int[1],
                                 conf.95.high = ft$conf.int[2],
                                 n = inft.inst,
                                 stringsAsFactors = FALSE)
      tmp2.annot.df <- rbind(tmp2.annot.df, tmp.annot.df)
    }
    bck.ovrl <- rowSums(as.data.frame(tmp.gr)[,setdiff(chrom.marks,mk)]) > 0
    for (chr.grp in chrom.marks) {
      if (!chr.grp %in% colnames(elementMetadata(tmp.gr))) {
        tmp.chr.df <- data.frame(chrom.state = chr.grp,
                                 peak.group = pk.grp, 
                                 mark = mk, 
                                 pvalue = -1,
                                 odds.ratio = NA,
                                 log2odds = NA,
                                 conf.95.low = NA,
                                 conf.95.high = NA,
                                 n = -1,
                                 fraction.total = 0,
                                 fraction.mark = 0,
                                 fraction.group = 0,
                                 peak.ratio = 0,
                                 stringsAsFactors = FALSE)
        tmp2.chr.df <- rbind(tmp2.chr.df, tmp.chr.df)  
        
      } else {
        inft.inst <- sum(elementMetadata(tmp.gr)[[chr.grp]] & tmp.gr$peak.group == pk.grp)
        
        inft.onst <- sum(elementMetadata(tmp.gr)[[chr.grp]] & (tmp.gr$peak.group != pk.grp))
        
        onft.inst <- sum((!elementMetadata(tmp.gr)[[chr.grp]]) & tmp.gr$peak.group == pk.grp)
        
        onft.onst <- sum(elementMetadata(tmp.gr)[[chr.grp]] | tmp.gr$peak.group == pk.grp)
        
        inft.bckg <- sum(elementMetadata(tmp.gr)[[chr.grp]] & bck.ovrl)
          
        am <- matrix(c(inft.inst,
                       inft.onst,
                       onft.inst,
                       length(tmp.gr) - onft.onst),
                     byrow = TRUE, ncol = 2)
        ft <- fisher.test(am)
        tmp.chr.df <- data.frame(chrom.state = chr.grp,
                                 peak.group = pk.grp, 
                                 mark = mk, 
                                 pvalue = ft$p.value,
                                 odds.ratio = ft$estimate,
                                 log2odds = log2(ft$estimate),
                                 conf.95.low = ft$conf.int[1],
                                 conf.95.high = ft$conf.int[2],
                                 n = inft.inst,
                                 fraction.total = inft.inst / length(tmp.gr),
                                 fraction.mark = inft.inst  / sum(elementMetadata(tmp.gr)[[chr.grp]]),
                                 fraction.group = inft.inst  / sum(tmp.gr$peak.group == pk.grp),
                                 peak.ratio = inft.inst / inft.bckg,
                                 stringsAsFactors = FALSE)
        tmp2.chr.df <- rbind(tmp2.chr.df, tmp.chr.df)  
      }
    }
  }
  
  tmp2.annot.df$padj <- p.adjust(tmp2.annot.df$pvalue)
  ft.annot.df <- rbind(ft.annot.df, tmp2.annot.df)
  
  tmp2.chr.df$padj <- p.adjust(tmp2.chr.df$pvalue)
  ft.chr.df <- rbind(ft.chr.df, tmp2.chr.df)
}

peak.set.lvls <- names(peak.set.cols.list$`7class`)

ft.chr.df$peak.group <- factor(ft.chr.df$peak.group,
                               levels = peak.set.lvls)

ft.chr.df$chrom.state <- factor(ft.chr.df$chrom.state,
                                levels = chrom.marks)
ft.chr.df$mark <- factor(ft.chr.df$mark, levels = c("H3K27me3","H3K9me3",
                                                    "H3K27ac","H3K4me3"))


p.cutoff <- 0.001

ft.chr.df$padj2 <- ft.chr.df$padj
ft.chr.df$log2odds2 <- ft.chr.df$log2odds
ft.chr.df$log2odds2[ft.chr.df$n <= 10 | ft.chr.df$padj2 >= p.cutoff] <- NA
ft.chr.df$label <- trimws(paste(ft.chr.df$n, pval2signf(ft.chr.df$padj2)))
ft.chr.df$peak.group2 <- as.character(ft.chr.df$peak.group)
ft.chr.df$peak.group2[ft.chr.df$peak.group == "Lost_3"] <- "Lost"
ft.chr.df$peak.group2[ft.chr.df$peak.group == "Gained_3"] <- "Gained"

ft.chr.df$peak.group2 <- factor(ft.chr.df$peak.group2, levels = c("Lost","Lost_2","Lost_1",
                                                                  "Kept",
                                                                  "Gained_1",
                                                                  "Gained_2","Gained"))
ft.chr.df$label[ft.chr.df$n < 10] <- as.character(ft.chr.df$n[ft.chr.df$n < 10])
ft.chr.df$label2 <- as.character(ft.chr.df$n)


ft.chr.df$percentage.total <- round(ft.chr.df$fraction.total * 100)
ft.chr.df$percentage.mark <- round(ft.chr.df$fraction.mark * 100)
ft.chr.df$percentage.group <- round(ft.chr.df$fraction.group * 100)

chip.Fig2C.plt <- filter(ft.chr.df, peak.group == "Lost_3") %>%
  ggplot(aes(y = chrom.state, x = mark, color = log2odds2, size = percentage.group)) +
  geom_point() + 
  scale_size_continuous(range = c(0, 12)) + 
  scale_color_distiller(type = "div", palette = "RdBu", limits = c(-2,2),oob=squish) +
  ylab("") + xlab("") +
  #guides(fill = guide_legend("Log-odds")) +
  labs(color = "Log-odds", size = "% Peaks in\n   Lost 3/3") +
  #geom_text(aes(label=label), size=4.5) +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  theme(panel.spacing = unit(0.3, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt))

ggsave(filename = file.path(plots.dir, "Fig2C_ChIP_lost_chrom_bubbleplot.pdf"),
       plot = chip.Fig2C.plt, width = 8, height = 10,
       device = cairo_pdf)


chip.FigS6A.plt <- filter(ft.chr.df, !chrom.state %in% c("uCpG.island")) %>%
  ggplot(aes(y = chrom.state, x = peak.group, 
             color = log2odds2, size = percentage.total)) +
  geom_point() + 
  scale_size_continuous(range = c(0, 12)) + 
  #scale_size_binned(range = c(0, 6)) +
  scale_color_distiller(type = "div", palette = "RdBu", limits = c(-2,2),oob=squish) +
  ylab("") + xlab("") +
  #guides(fill = guide_legend("Log-odds")) +
  labs(color = "Log-odds", size = "% of Query Peaks") +
  #geom_text(aes(label=label), size=4.5) +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  theme(panel.spacing = unit(0.3, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt)) +
  facet_wrap(~mark, nrow = 2)

ggsave(filename = file.path(plots.dir,"FigS6A_ChIP_peak_chromstate_bubbleplot.pdf"),
       plot = chip.FigS6A.plt, 
       device = cairo_pdf,
       width = 18, height = 12)

ft.annot.df$annotation <- factor(ft.annot.df$annotation,
                                 levels = names(annot.col))
ft.annot.df$padj2 <- ft.annot.df$padj
ft.annot.df$log2odds2 <- ft.annot.df$log2odds
ft.annot.df$log2odds2[ft.annot.df$n <= 10 | ft.annot.df$padj2 >= p.cutoff] <- NA

ft.annot.df$label <- trimws(paste(ft.annot.df$n, pval2signf(ft.annot.df$padj2)))
peak.set.lvls <- names(peak.set.cols.list$`7class`)
ft.annot.df$peak.group <- factor(ft.annot.df$peak.group,
                                 levels = peak.set.lvls)


chip.FigR8A.plt <- filter(ft.annot.df, mark == "H3K27me3") %>% 
  mutate(annotation = factor(annotation, 
                             levels = rev(c("Promoter","exon",
                                            "Distal Intergenic","intron",
                                            "5' UTR","3' UTR","Downstream <3kb")))) %>%
  ggplot(aes(y = annotation, x = peak.group, fill = log2odds2, label = label)) +
  geom_tile() + 
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-2,2)) +
  ylab("") + xlab("") +
  ggtitle("H3K27me3") +
  labs(fill = "Log-odds") +
  geom_text(aes(label=label), size=5) +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  theme(panel.spacing = unit(0.3, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt))

ggsave(filename = file.path(plots.dir,"FigR8A_ChIP_peaks_annot_H3K27me3.pdf"),
       plot = chip.FigR8A.plt, 
       device = cairo_pdf,
       width = 18, height = 12)




peak.widths.df <- c()

for (mk in mark.list) {
  merged.gr <- all.pk.mark.list[[mk]]
  pkgrps <- c("Lost_3","Kept","Gained_3")
  merged.gr$peak.group <- merged.gr$peak.7class

  lost.gr <- subset(merged.gr, peak.group == pkgrps[1])
  kept.gr <- subset(merged.gr, peak.group == "Kept")
  gained.gr <- subset(merged.gr, peak.group == pkgrps[3])
  
  for (cid in c("410","421","439","442")) {
    pkid <- paste0(mk, "_",cid)
    tmp.gr <- chip.peaks.list[[pkid]]
    tmp.gr <- tmp.gr[overlapsAny(tmp.gr, kept.gr, minoverlap = 1)]
    tmp.df <- data.frame(mark = mk, peak.group = "Kept",
                         clone = cid,
                         GC.content = tmp.gr$GC.cont,
                         width = width(tmp.gr), stringsAsFactors = FALSE)
    peak.widths.df <- rbind(peak.widths.df, tmp.df)
    
    
    if (cid != "410") {
      tmp.gr <- chip.peaks.list[[pkid]]
      tmp.gr <- tmp.gr[overlapsAny(tmp.gr, gained.gr, minoverlap = 1)]
      tmp.df <- data.frame(mark = mk, peak.group = "Gained",
                           clone = cid,
                           GC.content = tmp.gr$GC.cont,
                           width = width(tmp.gr), stringsAsFactors = FALSE)
      peak.widths.df <- rbind(peak.widths.df, tmp.df)
      
    } else {
      tmp.gr <- chip.peaks.list[[pkid]]
      tmp.gr <- tmp.gr[overlapsAny(tmp.gr, lost.gr, minoverlap = 1)]
      tmp.df <- data.frame(mark = mk, peak.group = "Lost",
                           clone = cid,
                           GC.content = tmp.gr$GC.cont,
                           width = width(tmp.gr), stringsAsFactors = FALSE)
      peak.widths.df <- rbind(peak.widths.df, tmp.df)
    }
  }
}

peak.widths.df$mark <- factor(peak.widths.df$mark)
peak.widths.df$clone <- factor(peak.widths.df$clone)
peak.widths.df$peak.group <- factor(peak.widths.df$peak.group, levels = c("Lost","Kept","Gained"))


chip.FigS5C.plt <- filter(peak.widths.df, !((peak.group == "Kept") & (clone != "410"))) %>%
  ggplot(aes(x = peak.group, y = log10(width), fill = peak.group)) +
  geom_violin() +
  geom_boxplot(fill = "darkgray", color = "black",
               width = 0.2) +
  theme_bw(base_size = 18) +
  ylab("log10(peak width)") +
  xlab("") +
  scale_fill_manual(values = peak.set.cols3) +
  facet_wrap(~mark, nrow = 2) +
  guides(fill=guide_legend(title="Peak group")) +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #legend.title = element_blank(),
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        strip.text.x = element_text(size = 6 * .pt),
        strip.text.y = element_text(size = 6 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3"),
        panel.spacing = unit(0.3, "lines"), 
        aspect.ratio = 0.8) #,

ggsave(filename = file.path(plots.dir, "FigS5C_ChIP_peak_widths.pdf"),
       device = cairo_pdf,plot = chip.FigS5C.plt,
       width = 14, height = 10)


### H3K27me3 peak overlaps with other sets
tmp.gr <- all.pk.mark.list$H3K27me3
tmp.gr$peak.group <- tmp.gr$peak.7class
tmp.gr$H2Aub <- factor(ifelse(tmp.gr$H2AK119ub1, "H2AK119ub1", "non-H2AK119ub1"),
                       levels = c("H2AK119ub1", "non-H2AK119ub1"))

h2aub.ft.df <- c()

for (pk.grp in unique(tmp.gr$peak.group)) {
  for (annot in unique(tmp.gr$H2Aub)) {
    inft.inst <- length(subset(tmp.gr, (H2Aub == annot) &
                                 (peak.group == pk.grp)))
    inft.onst <- length(subset(tmp.gr, (H2Aub == annot) &
                                 (peak.group != pk.grp)))
    onft.inst <- length(subset(tmp.gr, (H2Aub != annot) &
                                 (peak.group == pk.grp)))
    
    onft.onst <- length(subset(tmp.gr, (H2Aub == annot) |
                                 (peak.group == pk.grp)))
    
    am <- matrix(c(inft.inst,
                   inft.onst,
                   onft.inst,
                   length(tmp.gr) - onft.onst),
                 byrow = TRUE, ncol = 2)
    ft <- fisher.test(am)
    tmp.annot.df <- data.frame(annotation = annot,
                               peak.group = pk.grp, 
                               mark = mk, 
                               pvalue = ft$p.value,
                               odds.ratio = ft$estimate,
                               log2odds = log2(ft$estimate),
                               conf.95.low = ft$conf.int[1],
                               conf.95.high = ft$conf.int[2],
                               n = inft.inst,
                               stringsAsFactors = FALSE)
    h2aub.ft.df <- rbind(h2aub.ft.df, tmp.annot.df)
  }
}

h2aub.ft.df$padj <- p.adjust(h2aub.ft.df$pvalue, method = "BH")

h2aub.ft.df$label <- paste(h2aub.ft.df$n, pval2signf(h2aub.ft.df$padj))
h2aub.ft.tmp.df <- filter(h2aub.ft.df, peak.group %in% c("Gained_3","Kept","Lost_3"))

h2aub.ft.tmp.df$peak.group[h2aub.ft.tmp.df$peak.group == pkgrps[1]] <- "Lost"
h2aub.ft.tmp.df$peak.group[h2aub.ft.tmp.df$peak.group == pkgrps[3]] <- "Gained"

h2aub.ft.tmp.df$peak.group <- factor(h2aub.ft.tmp.df$peak.group, levels = c("Lost","Kept", "Gained"))
h2aub.ft.tmp.df$annotation <- factor(h2aub.ft.tmp.df$annotation, levels = c("H2AK119ub1", "non-H2AK119ub1"))
h2aub.ft.tmp.df$log2odds[h2aub.ft.tmp.df$padj > 0.001] <- NA

chip.figS6B.plt <- ggplot(h2aub.ft.tmp.df, aes(x = peak.group, y = annotation, 
                                               fill = log2odds, label = label)) +
  geom_tile() + 
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-2,2), oob=squish) +
  ylab("") + xlab("H3K27me3 peaks") +
  labs(fill = "Log-odds") +
  geom_text(aes(label=label), size=8) +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  theme(panel.spacing = unit(0.3, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt))

ggsave(filename = file.path(plots.dir, "FigS6B_ChIP_K27me3_H2Aub_peak_overlaps_heatmap.pdf"),
       device = cairo_pdf,plot = chip.figS6B.plt,
       width = 12, height = 8)






### SUZ12 analysis
percentageOverlapping <- function(gr1, gr2, percent = 0.5) {
  #gr1 <- chip.cond.lenient.overlap.peaks.list$H3K4me3$WT
  #gr2 <- chip.cond.lenient.overlap.peaks.list$H3K27me3$WT
  
  hits <- findOverlaps(gr1, gr2)
  overlaps <- pintersect(gr1[queryHits(hits)], gr2[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(gr1[queryHits(hits)])
  hits <- hits[percentOverlap >= percent]
  return(length(unique(subjectHits(hits))) / length(gr1))
}

isOverlapping <- function(gr1, gr2, percent = 0.5) {
  hits <- findOverlaps(gr1, gr2)
  overlaps <- pintersect(gr1[queryHits(hits)], gr2[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(gr1[queryHits(hits)])
  hits <- hits[percentOverlap >= percent]
  return(names(gr1) %in% names(gr1[queryHits(hits)]))
}

percentageCovering <- function(gr1, gr2, percent = 0.5) {
  hits <- findOverlaps(gr1, gr2)
  overlaps <- pintersect(gr1[queryHits(hits)], gr2[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(gr2[subjectHits(hits)])
  hits <- hits[percentOverlap >= percent]
  return(length(unique(subjectHits(hits))) / length(gr1))
}

isCovering <- function(gr1, gr2, percent = 0.5) {
  hits <- findOverlaps(gr1, gr2)
  overlaps <- pintersect(gr1[queryHits(hits)], gr2[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(gr2[subjectHits(hits)])
  hits <- hits[percentOverlap >= percent]
  return(names(gr1) %in% names(gr1[queryHits(hits)]))
}


### SUZ12 Nucleation sites
suz12.strict.peak.list <- chip.peaks.list[paste0("SUZ12_",c("410","439","442"))]
suz12.strict.peak.list <- lapply(suz12.strict.peak.list, FUN = function(x) {
  return(subset(x, qval > -log10(0.05)))
})

wt.k27.gr <- chip.peaks.list[["H3K27me3_410"]]
wt.suz12.gr <- suz12.strict.peak.list[["SUZ12_410"]]
suz12.K27me3.ovrl.df <- c()
for (cid in c("439","442")) {
  mut.k27.gr <- chip.peaks.list[[paste0("H3K27me3_",cid)]]
  mut.suz12.gr <- suz12.strict.peak.list[[paste0("SUZ12_",cid)]]
  
  k27.ovrl <- findOverlapsOfPeaks(list("Peaks_WT" = wt.k27.gr, 
                                       "Peaks_MUT" = mut.k27.gr),
                                  minoverlap = 1)
  suz12.ovrl <- findOverlapsOfPeaks(list("Peaks_WT" = wt.suz12.gr, 
                                         "Peaks_MUT" = mut.suz12.gr),
                                  minoverlap = 1)
  
  k27.lost.gr <- k27.ovrl$peaklist$Peaks_WT
  k27.lost.gr$peak.group <- "Lost"
  k27.kept.gr <- k27.ovrl$peaklist$`Peaks_WT///Peaks_MUT`
  k27.kept.gr$peak.group <- "Kept"
  k27.gained.gr <- k27.ovrl$peaklist$Peaks_MUT
  k27.gained.gr$peak.group <- "Gained"
  all.k27.gr <- c(k27.lost.gr, c(k27.kept.gr, k27.gained.gr))
  
  suz12.lost.gr <- suz12.ovrl$peaklist$Peaks_WT
  suz12.lost.gr$peak.group <- "Lost"
  suz12.kept.gr <- suz12.ovrl$peaklist$`Peaks_WT///Peaks_MUT`
  suz12.kept.gr$peak.group <- "Kept"
  suz12.gained.gr <- suz12.ovrl$peaklist$Peaks_MUT
  suz12.gained.gr$peak.group <- "Gained"
  all.suz12.gr <- c(suz12.lost.gr, c(suz12.kept.gr, suz12.gained.gr))
  #universe.gr <- all.suz12.gr
  universe.gr <- GenomicRanges::reduce(c(chip.peaks.list[["SUZ12_410"]], 
                                         chip.peaks.list[[paste0("SUZ12_",cid)]]))
  
  for (dt.suz12 in c("Lost","Kept","Gained")) {
    suz12.gr <- subset(all.suz12.gr, peak.group == dt.suz12)
    names(suz12.gr) <- paste0("SUZ12_", 1:length(suz12.gr))
    for (dt.k27me3 in c("Lost","Kept","Gained")) {
      k27me3.gr <- subset(all.k27.gr, peak.group == dt.k27me3)
      names(k27me3.gr) <- paste0("SUZ12_", 1:length(k27me3.gr))
      ovrl <- findOverlapsOfPeaks(list("SUZ12" = suz12.gr, "H3K27me3" = k27me3.gr),
                                  minoverlap = 1, connectedPeaks = "merge")
      #suz12.and.k27me3 <- suz12.gr[isOverlapping(suz12.gr,k27me3.gr, percent = 0.8)]
      #suz12.not.k27me3 <- suz12.gr[!isOverlapping(suz12.gr,k27me3.gr, percent = 0.8)]
      #k27me3.not.suz12 <- k27me3.gr[!isCovering(k27me3.gr, suz12.gr, percent = 0.8)]
      suz12.and.k27me3 <- ovrl$peaklist$`SUZ12///H3K27me3`
      suz12.not.k27me3 <- ovrl$peaklist$SUZ12
      k27me3.not.suz12 <- ovrl$peaklist$H3K27me3
      
      the.rest <- universe.gr[!(overlapsAny(universe.gr,suz12.gr, minoverlap = 1) | 
                                  overlapsAny(universe.gr,k27me3.gr, minoverlap = 1))]
      tb <- matrix(c(length(suz12.and.k27me3), length(k27me3.not.suz12), 
                     length(suz12.not.k27me3), length(the.rest)), nrow = 2)
      ft <- fisher.test(tb)
      tmp.df <- data.frame(clone = cid,
                           SUZ12.group = dt.suz12, 
                           H3K27me3.group = dt.k27me3, 
                           pvalue = ft$p.value, odds = ft$estimate,
                           lower.conf = ft$conf.int[1], upper.conf = ft$conf.int[2])
      tmp.df$n <- length(suz12.and.k27me3)
      suz12.K27me3.ovrl.df <- rbind(suz12.K27me3.ovrl.df, tmp.df)
    }
  }
}

mx.odds <- min(max(abs(log2(suz12.K27me3.ovrl.df$odds))),2)

suz12.K27me3.ovrl.df$padj <- p.adjust(suz12.K27me3.ovrl.df$pvalue, method = "BH")
suz12.K27me3.ovrl.df$label <- trimws(paste(suz12.K27me3.ovrl.df$n, pval2signf(suz12.K27me3.ovrl.df$padj)))
suz12.K27me3.ovrl.df$SUZ12.group <- factor(suz12.K27me3.ovrl.df$SUZ12.group, levels = c("Lost","Kept","Gained"))
suz12.K27me3.ovrl.df$H3K27me3.group <- factor(suz12.K27me3.ovrl.df$H3K27me3.group, levels = c("Lost","Kept","Gained"))
suz12.K27me3.ovrl.df$clone <- factor(suz12.K27me3.ovrl.df$clone)

suz12.K27me3.ovrl.df$logodds <- log2(suz12.K27me3.ovrl.df$odds)
suz12.K27me3.ovrl.df$logodds[suz12.K27me3.ovrl.df$padj > 0.001] <- NA

chip.FigS8C.plt <- suz12.K27me3.ovrl.df %>% mutate(logodds = log(odds)) %>%
  mutate(logodds = ifelse(logodds > 2, 2, logodds)) %>%
  mutate(logodds = ifelse(logodds < -2, -2, logodds)) %>%
  mutate(logodds = ifelse(padj > p.cutoff, NA, logodds)) %>%
  ggplot(aes(x = SUZ12.group, y = H3K27me3.group, fill = logodds, label = label)) +
  geom_tile() +
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-mx.odds,mx.odds)) +
  ylab("H3K27me3 peaks") + xlab("SUZ12 peaks") +
  labs(fill = "Log-odds") +
  geom_text(aes(label=label), size=7) +
  theme_bw(base_family = "Helvetica", base_size = 18) +
  theme(panel.spacing = unit(0.3, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        aspect.ratio = 0.8) +
  facet_wrap(~clone, labeller = labeller(.cols = clone.labs))

ggsave(filename = file.path(plots.dir,"FigS8C_ChIP_K27me3_SUZ12_peak_overlaps_heatmap.pdf"),
       device = cairo_pdf,plot = chip.FigS8C.plt,
       width = 12, height = 10)




## Distance based analysis between K27me3 and SUZ12

all.k27me3.gr <- all.pk.mark.list$H3K27me3
all.k27me3.gr$peak.group <- all.k27me3.gr$peak.7class

k27.lost.gr <- subset(all.k27me3.gr, peak.group == pkgrps[1])
k27.lost.gr$peak.set <- "Lost"
k27.kept.gr <- subset(all.k27me3.gr, peak.group == "Kept")
k27.kept.gr$peak.set <- "Kept"
k27.gained.gr <- subset(all.k27me3.gr, peak.group == pkgrps[3])
k27.gained.gr$peak.set <- "Gained"

tmp.k27.gr <- c(k27.lost.gr, c(k27.kept.gr, k27.gained.gr))


wt.suz12.gr <- suz12.strict.peak.list$SUZ12_410
names(wt.suz12.gr) <- paste0("Peak_", 1:length(wt.suz12.gr))


wt.suz12.gr$H3K27me3.group <- factor(tmp.k27.gr$peak.set[nearest(resize(wt.suz12.gr, width = 1, fix = "center"), 
                                                                 resize(tmp.k27.gr, width = 1, fix = "center"))], levels = c("Lost","Kept","Gained"))

wt.suz12.gr <- subset(wt.suz12.gr, !is.na(H3K27me3.group))

tmp.k27.lost.gr <- subset(all.k27me3.gr, peak.group == pkgrps[1])
tmp.k27.lost.gr <- resize(tmp.k27.lost.gr, fix = "center", width = 1)
tmp.k27.kept.gr <- subset(all.k27me3.gr, peak.group == "Kept")
tmp.k27.kept.gr <- resize(tmp.k27.kept.gr, fix = "center", width = 1)

tmp.suz12.gr <- suz12.strict.peak.list$SUZ12_410
tmp.suz12.gr <- resize(tmp.suz12.gr, width = 1, fix = "center")

dist.suz12.lost <- distanceToNearest(tmp.k27.lost.gr, tmp.suz12.gr)
dist.suz12.kept <- distanceToNearest(tmp.k27.kept.gr, tmp.suz12.gr)

k27.suz12.dist.all.df <- rbind(data.frame(peak.group = "Lost", 
                                          kb = (start(tmp.k27.lost.gr[queryHits(dist.suz12.lost)]) -
                                                  start(tmp.suz12.gr[subjectHits(dist.suz12.lost)])) / 1000,
                                          stringsAsFactors = FALSE),
                               data.frame(peak.group = "Kept", 
                                          kb = (start(tmp.k27.kept.gr[queryHits(dist.suz12.kept)]) -
                                                  start(tmp.suz12.gr[subjectHits(dist.suz12.kept)])) / 1000,
                                          stringsAsFactors = FALSE))

k27.suz12.dist.all.df$peak.group <- factor(k27.suz12.dist.all.df$peak.group, levels = c("Lost","Kept"))


all.k9me3.gr <- all.pk.mark.list$H3K9me3
all.k9me3.gr$peak.group <- all.k9me3.gr$peak.7class

k9.lost.gr <- subset(all.k9me3.gr, peak.group == pkgrps[1])
k9.lost.gr$peak.set <- "Lost"
k9.kept.gr <- subset(all.k9me3.gr, peak.group == "Kept")
k9.kept.gr$peak.set <- "Kept"
k9.gained.gr <- subset(all.k9me3.gr, peak.group == pkgrps[3])
k9.gained.gr$peak.set <- "Gained"

tmp.k9.gr <- c(k9.lost.gr, c(k9.kept.gr, k9.gained.gr))


k9f.top.list <- list()
k9f.dist.all.df <- c()
k9f.density.df <- c()
for (k9f in c("SETDB1","SUV39H1","SUV39H2")) {
  top.gr <- chip.external.peaks.list[[k9f]]
  top.gr <- top.gr[order(top.gr$pval, top.gr$log2FC, decreasing = TRUE)][1:10000]
  k9f.top.list[[k9f]] <- top.gr
  top.gr <- resize(top.gr, width = 1, fix = "center")
  
  tmp.lost.gr <- subset(all.k9me3.gr, peak.group == pkgrps[1])
  tmp.lost.gr <- resize(tmp.lost.gr, width = 1, fix = "center")
  tmp.kept.gr <- subset(all.k9me3.gr, peak.group == "Kept")
  tmp.kept.gr <- resize(tmp.kept.gr, width = 1, fix = "center")
  
  dist.lost <- distanceToNearest(tmp.lost.gr, top.gr)
  dist.kept <- distanceToNearest(tmp.kept.gr, top.gr)
  
  tmp.df <- rbind(data.frame(peak.group = "Lost", 
                             kb = (start(tmp.lost.gr[queryHits(dist.lost)]) -
                                     start(top.gr[subjectHits(dist.lost)])) / 1000,
                             stringsAsFactors = FALSE),
                  data.frame(peak.group = "Kept", 
                             kb = (start(tmp.kept.gr[queryHits(dist.kept)]) -
                                     start(top.gr[subjectHits(dist.kept)])) / 1000,
                             stringsAsFactors = FALSE))
  tmp.df$K9.factor <- k9f
  k9f.dist.all.df <- rbind(k9f.dist.all.df, tmp.df)
  tlost <- density(filter(tmp.df, peak.group == "Lost")$kb)
  densisty.df <- data.frame(dx = tlost$x, dy = tlost$y, 
                            peak.group = "Lost", K9.factor = k9f, stringsAsFactors = FALSE)
  k9f.density.df <- rbind(k9f.density.df, densisty.df)
  tkept <- density(filter(tmp.df, peak.group == "Kept")$kb)
  densisty.df <- data.frame(dx = tkept$x, dy = tkept$y, 
                            peak.group = "Kept", K9.factor = k9f, 
                            stringsAsFactors = FALSE)
  k9f.density.df <- rbind(k9f.density.df, densisty.df)
  
}

k9f.dist.all.df$peak.group <- factor(k9f.dist.all.df$peak.group, levels = c("Lost","Kept"))

suz12.tmp.ovrl <- k27.suz12.dist.all.df
suz12.tmp.ovrl$factor <- "SUZ12"

setdb.tmp.ovrl <- k9f.dist.all.df
setdb.tmp.ovrl$factor <- as.character(k9f.dist.all.df$K9.factor)
setdb.tmp.ovrl <- setdb.tmp.ovrl[,c("peak.group","kb","factor")]

all.dist.fact.df <- rbind(suz12.tmp.ovrl, setdb.tmp.ovrl)
all.dist.fact.df$factor <- factor(all.dist.fact.df$factor,
                                  levels = c("SUZ12","SETDB1","SUV39H1","SUV39H2"))



chip.Fig2C.plt <- filter(all.dist.fact.df, factor != "SETDB1") %>%
  ggplot(aes(x=abs(kb), fill=peak.group, color = peak.group, alpha = peak.group)) +
  geom_density(size = 1.3) +
  scale_fill_manual(values = c("Lost" = "darkblue", "Kept" = "darkgray")) +
  scale_color_manual(values = c("Lost" = "darkblue", "Kept" = "darkgray")) +
  scale_alpha_manual(values = c("Lost" = 0.5, "Kept" = 0.3)) +
  guides(fill =  guide_legend(title = "peaks"),
         color =  guide_legend(title = "peaks"),
         alpha =  "none") +
  xlim(0,500) +
  ylab("Density of peaks") +
  xlab("Distance (kb) to nearest peak center ") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 5 * .pt),
        legend.text = element_text(size = 6 * .pt),
        legend.title = element_text(size = 6 * .pt),
        panel.spacing = unit(0.35, "lines"),
        strip.text.x = element_text(size = 5 * .pt),
        strip.text.y = element_text(size = 5 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3")) +
  facet_wrap(~factor, nrow = 1, scales = "free_y")


ggsave(filename = file.path(plots.dir,"Fig2C_ChIP_fact_distance_density.pdf"),
       plot = chip.Fig2C.plt,  
       width = 12,  height = 4)



chip.FigS5B.plt <- filter(all.dist.fact.df, factor == "SETDB1") %>%
  ggplot(aes(x=abs(kb), fill=peak.group, color = peak.group, alpha = peak.group)) +
  geom_density(size = 1.3) +
  scale_fill_manual(values = c("Lost" = "darkblue", "Kept" = "darkgray")) +
  scale_color_manual(values = c("Lost" = "darkblue", "Kept" = "darkgray")) +
  scale_alpha_manual(values = c("Lost" = 0.5, "Kept" = 0.3)) +
  guides(fill =  guide_legend(title = "H3K9me3\n peaks"),
         color =  guide_legend(title = "H3K9me3\n peaks"),
         alpha =  "none") +
  xlim(0,500) +
  ylab("Density of peaks") +
  xlab("Distance (kb) to nearest SETDB1 peak") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 5 * .pt),
        legend.text = element_text(size = 6 * .pt),
        legend.title = element_text(size = 6 * .pt),
        panel.spacing = unit(0.35, "lines"),
        strip.text.x = element_text(size = 5 * .pt),
        strip.text.y = element_text(size = 5 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3"))


ggsave(filename = file.path(plots.dir,"FigS5B_ChIP_SETDB1_H3K9me3_density.pdf"),
       plot = chip.FigS5B.plt,  
       width = 14,  height = 12)


suz12.olaps <- ssvOverlapIntervalSets(suz12.strict.peak.list, minoverlap = 1)

chip.figS8B.plt <- ssvFeatureEuler(suz12.olaps, circle_colors = c("blue","orange","red")) + #, counts_txt_size = 6) + 
  theme_void(base_size = 20, base_family = "Helvetica") + 
  theme(legend.position = "top", aspect.ratio = 0.9)

ggsave(filename = file.path(plots.dir,"FigS8B_SUZ12_peak_euler.pdf"),
       plot = chip.figS8B.plt, 
       width = 12,  height = 10)



### ==== qChIP analyses =======================================================
read.param <- readParam(discard = mm10.blacklist.noovrl.gr,
                        restrict = chroms, pe = "first")

chip.bam.dir <- "data/ChIP_seq/MAPPING_mm10/bam_files"
chip.ext.bam.dir <- "data/ChIP_seq/external/MAPPING_mm10/bam_files"

cl <- makeCluster(min(5, num.cores), 
                  type = "FORK")


registerDoParallel(cl)

all.csaw.chip.windows.rrpm.mean.df <- foreach(mk = c(mark.list,"SUZ12"),
                                              .combine = "rbind",
        .packages = c("tidyverse","csaw", "edgeR","limma")) %dopar% {
  #print(mk)
  ws <- 1000
  if (mk %in% c("H3K9me3","H3K27me3")) {
    ws <- 5000
  }
  load(file.path(chip.data.dir, paste0("ChIP_csaw_windows_", mk,"_w",ws,".RData")))
  tmp.stats.df <- filter(chip.stats.df, mark == mk)
  tmp.stats.df <- tmp.stats.df[colnames(csaw.chip.windows),]
  
  csaw.windows.log2FC.mat <- csaw.windows.log2FC.mat[rownames(csaw.chip.windows),]
  csaw.windows.filter.mat <- csaw.windows.log2FC.mat >= log2(1.5)
  
  csaw.chip.windows <- csaw.chip.windows[rowSums(csaw.windows.filter.mat) > 0,]
  csaw.chip.windows.rrpm.df <- cpm(assay(csaw.chip.windows), 
                                   lib.size = csaw.chip.windows$exo.reads,
                                   log = TRUE)
  
  csaw.chip.windows.rrpm.melt.df <- melt(csaw.chip.windows.rrpm.df)
  colnames(csaw.chip.windows.rrpm.melt.df) <- c("name","ID","logRRPM")
  csaw.chip.windows.rrpm.melt.df$clone <- strsplit2(as.character(csaw.chip.windows.rrpm.melt.df$ID),
                                                    "_")[,1]

  csaw.chip.windows.rrpm.mean.df <- csaw.chip.windows.rrpm.melt.df %>% group_by(name,clone) %>%
    summarise(logRRPM = mean(logRRPM, na.rm=T)) %>% 
    mutate(condition = factor(gsub("-","_", clone.conds[as.character(clone)]),
                              levels = c("WT","MCM2_2A")),
           mark = mk) %>% as.data.frame()
  return(csaw.chip.windows.rrpm.mean.df)
  
}

stopCluster(cl)


tmp.df <- filter(all.csaw.chip.windows.rrpm.mean.df, mark != "SUZ12",
                 clone %in% c("410","421","439","442")) %>%
  group_by(mark, clone) %>% mutate(Q3 = quantile(logRRPM, 0.75), 
                                   Q1 = quantile(logRRPM, 0.25), 
                                   IQR = IQR(logRRPM)) %>%
  filter(logRRPM < (Q3 + 1.5*IQR),  logRRPM > (Q1 - 1.5 * IQR)) %>%
  mutate(min.RFD = min(logRRPM, na.rm = TRUE),
         max.RFD = max(logRRPM, na.rm = TRUE)) %>% select(clone,mark, 
                                                          min.RFD, max.RFD) %>%
  unique() %>% as.data.frame()

tmp.df$mark <- factor(tmp.df$mark, levels = mark.list)

mark.minRRPM <- lapply(levels(tmp.df$mark), FUN = function(mk) {
  return(min(filter(tmp.df, mark == mk)$min.RFD) - 1)  
})

names(mark.minRRPM) <- levels(tmp.df$mark)

mark.maxRRPM <- lapply(levels(tmp.df$mark), FUN = function(mk) {
  return(max(filter(tmp.df, mark == mk)$max.RFD) + 1)  
})

names(mark.maxRRPM) <- levels(tmp.df$mark)


chip.FigS7h.plt <- 
  filter(all.csaw.chip.windows.rrpm.mean.df, mark != "SUZ12",
         clone %in% c("410","439","442","421")) %>%
  mutate(mark = factor(mark, levels = c("H3K27me3","H3K9me3","H3K27ac","H3K4me3")),
         clone = factor(clone, levels = c("410","439","442","421"))) %>%
  ggplot(aes(x = clone, y = logRRPM, fill = condition)) +
  geom_boxplot(outlier.shape = NA,
               notch = 1, lwd = 0.8, width = 0.5) + 
  scale_fill_manual(values=c("WT" = "blue", "MCM2_2A" = "orange")) +
  ylab("qChIP signal [log10(RRPM)]") +
  xlab("Cell line") +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  scale_x_discrete(labels = clone.labs[c("410","421","439","442")]) +
  guides(fill = guide_legend(title = "")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        axis.text.x = element_text(size = 6 * .pt, angle = 45, hjust = 1),
        legend.text = element_text(size = 6 * .pt),
        panel.spacing = unit(0.35, "lines"),
        strip.text.x = element_text(size = 6 * .pt),
        strip.text.y = element_text(size = 6 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3")) +
  facet_wrap(~mark, nrow = 2, scales = "free_y") +
  facetted_pos_scales(y = list(mark == "H3K27ac" ~ scale_y_continuous(limits = c(mark.minRRPM[["H3K27ac"]],
                                                                                 mark.maxRRPM[["H3K27ac"]])),
                               mark == "H3K27me3" ~ scale_y_continuous(limits = c(mark.minRRPM[["H3K27me3"]],
                                                                                  mark.maxRRPM[["H3K27me3"]])),
                               mark == "H3K4me3" ~ scale_y_continuous(limits = c(mark.minRRPM[["H3K4me3"]],
                                                                                 mark.maxRRPM[["H3K4me3"]])),
                               mark == "H3K9me3" ~ scale_y_continuous(limits = c(mark.minRRPM[["H3K9me3"]],
                                                                                 mark.maxRRPM[["H3K9me3"]]))))


ggsave(filename = file.path(plots.dir,"FigS7H_ChIP_marks_levels_boxplots.pdf"),
       plot = chip.FigS7h.plt, 
       device = cairo_pdf,
       width = 12,  height = 10)




tmp2.df <- filter(all.csaw.chip.windows.rrpm.mean.df, mark == "SUZ12",
                 clone %in% c("410","439","442")) %>%
  group_by(clone) %>% mutate(Q3 = quantile(logRRPM, 0.75), 
                                   Q1 = quantile(logRRPM, 0.25), 
                                   IQR = IQR(logRRPM)) %>%
  filter(logRRPM < (Q3 + 1.5*IQR),  logRRPM > (Q1 - 1.5 * IQR)) %>%
  mutate(min.RFD = min(logRRPM, na.rm = TRUE),
         max.RFD = max(logRRPM, na.rm = TRUE)) %>% select(clone,
                                                          min.RFD, max.RFD) %>%
  unique() %>% as.data.frame()




chip.FigS8a.plt <- 
  filter(all.csaw.chip.windows.rrpm.mean.df, mark == "SUZ12",
         clone %in% c("410","439","442")) %>%
  mutate(clone = factor(clone, levels = c("410","439","442"))) %>%
  ggplot(aes(x = clone, y = logRRPM, fill = condition)) +
  geom_boxplot(outlier.shape = NA,
               notch = 1, lwd = 0.8, width = 0.5) + 
  scale_fill_manual(values=c("WT" = "blue", "MCM2_2A" = "orange")) +
  ylab("SUZ12 qChIP signal [log10(RRPM)]") +
  xlab("Cell line") +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  scale_x_discrete(labels = clone.labs[c("410","439","442")]) +
  guides(fill = guide_legend(title = "")) +
  ylim(min(tmp2.df$min.RFD) - 0.5,
       max(tmp2.df$max.RFD) + 0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        axis.text.x = element_text(size = 6 * .pt, angle = 45, hjust = 1),
        legend.text = element_text(size = 6 * .pt),
        panel.spacing = unit(0.35, "lines"),
        strip.text.x = element_text(size = 6 * .pt),
        strip.text.y = element_text(size = 6 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3")) 

ggsave(filename = file.path(plots.dir,"FigS8A_ChIP_SUZ12_levels_boxplots.pdf"),
       plot = chip.FigS7h.plt, 
       device = cairo_pdf,
       width = 12,  height = 10)



### Correlation between ChIP rescue samples
load("data/ChIP_seq/ChIP_rescue_csaw_windows_H3K27me3_w1000.RData")


chip.rescue.ovl <- findOverlapsOfPeaks(chip.rescue.peaks.list,
                                       minoverlap = 1, 
                                       connectedPeaks = "merge")

chip.rescue.cons.gr <- GRanges()
pkl <- names(chip.rescue.ovl$peaklist)
for (j in pkl) {
  chip.rescue.cons.gr <- c(chip.rescue.cons.gr,chip.rescue.ovl$peaklist[[j]])
}

keep.windows <- overlapsAny(rowRanges(csaw.chip.windows),
                            chip.rescue.cons.gr, type = "within")

csaw.dge <- asDGEList(csaw.chip.windows[keep.windows,])
csaw.dge <- calcNormFactors(csaw.dge, method = "TMM")
csaw.dge.tmm <- cpm(csaw.dge)
colnames(csaw.dge.tmm) <- colnames(csaw.chip.windows)


chip.rg.cor.df <- cor(csaw.dge.tmm, method = "spearman")

tmp.annot.df <- chip.rescue.stats.df[rownames(chip.rg.cor.df),]

rownames(tmp.annot.df) <- paste(clone.labs[as.character(tmp.annot.df$clone)],tmp.annot.df$replicate)
colnames(chip.rg.cor.df) <- rownames(tmp.annot.df)
rownames(chip.rg.cor.df) <- rownames(tmp.annot.df)


annot.cols <- list(condition = c("WT" = "darkblue", 
                                 "MCM2-2A" = "orange",
                                 "MCM2-R" = "purple"))

colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

cairo_pdf(filename = file.path(plots.dir, "FigS8F_ChIP_rescue_K27me3_cor_heatmap.pdf"),
          width = 18, height = 16)

phout <- pheatmap(chip.rg.cor.df,
                  col = colors,
                  fontsize_col = 14,
                  fontsize_row = 14,
                  annotation_colors = annot.cols,
                  annotation_row = tmp.annot.df[,"condition",drop=F],
                  annotation_col = tmp.annot.df[,"condition",drop=F],
                  cellwidth = 40, cellheight = 40,
                  fontsize = 12)

print(phout)
dev.off()



#### Stats for qChIP

qchip.mark.list <- lapply(unique(all.csaw.chip.windows.rrpm.mean.df$mark), FUN = function(mk) {
  return(filter(all.csaw.chip.windows.rrpm.mean.df, mark == mk))
})

names(qchip.mark.list) <- unique(all.csaw.chip.windows.rrpm.mean.df$mark)


wilcox.comparisons <- list(c("410","439"),c("410","442"),c("410","421"),
                           c("439","442"),c("439","421"),c("442","421"))


mk.df <- qchip.mark.list$H3K27me3 %>%
  mutate(clone = factor(clone, levels = c("410","439","442","421")))

wilcox.es.test.df <- wilcox_effsize(formula = logRRPM ~ clone,
                                       data = qchip.mark.list$H3K27me3 %>%
                                         mutate(clone = factor(clone, levels = c("410","439","442","421"))),
                                      comparisons = wilcox.comparisons,
                                 paired = FALSE) %>% as.data.frame()


tmp.df <- compare_means(formula = logRRPM ~ clone,
                        data = mk.df,
                        method = "wilcox.test")


wilcox.test.df <- full_join(compare_means(formula = logRRPM ~ clone,
                                          data = mk.df,
                                          method = "wilcox.test") %>%
                              as.data.frame(), 
                            wilcox_effsize(logRRPM ~ clone, data = mk.df) %>%
                            as.data.frame())

cl <- makeCluster(5, type = "FORK")
registerDoParallel(cl)
all.qchip.stats.df <-  foreach(df = iter(qchip.mark.list), .combine = "rbind", 
                               .packages = c("coin","rstatix", "tidyverse",
                                             "ggpubr")) %dopar% {
                                 mk <- unique(df$mark)
                                 wilcox.comparisons <- list(c("410","439"),c("410","442"),c("410","421"),
                                                            c("439","442"),c("439","421"),c("442","421"))
                                 lvls <- c("410","439","442","421")
                                 if (mk == "SUZ12") {
                                   wilcox.comparisons <- list(c("410","439"),c("410","442"),
                                                              c("439","442"))
                                   lvls <- c("410","439","442")
                                 }
                                 df$clone <- factor(df$clone, levels = lvls)
                                 
                                 wilcox.test.df <- full_join(compare_means(formula = logRRPM ~ clone,
                                                                           data = df,
                                                                           method = "wilcox.test") %>%
                                                               as.data.frame(), 
                                                             wilcox_effsize(logRRPM ~ clone, data = df) %>%
                                                               as.data.frame()) %>%
                                   mutate(mark = mk)
                                 
                                 return(wilcox.test.df)
                               }

stopCluster(cl)

write.xlsx(all.qchip.stats.df,
           file = file.path(chip.data.dir, 
                            paste0("qChIP_stats.xlsx")))
