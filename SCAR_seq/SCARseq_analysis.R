library(GenomicAlignments)
library(GenomicFeatures)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(ChIPseeker)
library(stringi)
library(clusterProfiler)
library(foreach)
library(doParallel)
library(BiocParallel)
library(ggplot2)
library(ggthemr)
library(ggrepel)
library(ggpubr)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(gridExtra)
library(scales)
library(pheatmap)
library(reshape2)
library(bamsignals)
library(ggpmisc)
library(csaw)
library(patchwork)
library(extrafont)
library(remotes)
library(ecoflux)
library(factoextra)
library(limma)
library(tidyverse)
library(latex2exp)
library(ggsci)
library(edgeR)
library(coin)
library(rstatix)
library(openxlsx)


### Set working directory and parameters ###


wd.dir <- "Wenger_et_al_2022"
setwd(wd.dir)

data.dir <- "data"
plots.dir <- "plots"

chip.data.dir <- file.path(data.dir,"ChIP_seq")
scar.data.dir <- file.path(data.dir,"SCAR_seq")

external.data.dir <- "data/external"


chroms <- paste("chr",1:19,sep="") ## For now we don't work with sex chromosomes
spike.chroms <- c("chr2L_spike", "chr2LHet_spike","chr2R_spike", "chr2RHet_spike", 
                  "chr3L_spike", "chr3LHet_spike", "chr3R_spike", "chr3RHet_spike", 
                  "chr4_spike", "chrM_spike", "chrU_spike", "chrUextra_spike",
                  "chrX_spike", "chrXHet_spike","chrYHet_spike")

## ===================  Load external data which will be used later ============
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



#### Define colors
mk.tp.pal.list <- list()

mk.tp.pal.list[["WT"]] <- c(pal_material("blue")(10)[c(10,7,5,3)],"#043927")
mk.tp.pal.list[["MCM2-2A"]] <- c(pal_material("orange")(10)[c(10,7,5,3)],"#043927")

names(mk.tp.pal.list[["MCM2-2A"]]) <- c("T=0","T=1","T=3","T=8","OK-seq")
names(mk.tp.pal.list[["MCM2-2A"]]) <- c("T=0","T=1","T=3","T=8","OK-seq")


part.cond.cols <- c("OK-seq" = "#043927",
                    "WT" = c(pal_material("blue")(10)[10]),
                    "MCM2-2A" = c(pal_material("orange")(10)[5]),
                    "MCM2-Y90A" = c(pal_material("orange")(10)[3]),
                    "MCM2-2A-R" = c(pal_material("purple")(10)[10]),
                    "POLE4-KO" = c(pal_material("red")(10)[10])) 

tp.linetype.list <- c("T=0" = "solid", "T=1" = "dotdash", "T=3" = "dashed", "T=8" = "dotted")
tp.alpha.list <- c("T=0" = 1, "T=1" = 0.7, "T=3" = 0.4, "T=8" = 0.2)

tp.pal.list <- list()
tp.pal.list[["H3K27me3"]] <- brewer.pal(n=9,name="Purples")[c(8,6,4)]
tp.pal.list[["H3K4me3"]] <- brewer.pal(n=9,name="Greens")[c(8,6,4)]
tp.pal.list[["H3K9me3"]] <- brewer.pal(n=9,name="Oranges")[c(8,6,4)]
tp.pal.list[["H3K27ac"]] <- brewer.pal(n=9,name="Blues")[c(9,7,5,3)]

names(tp.pal.list[["H3K27me3"]]) <- c("T=0", "T=3", "T=8")
names(tp.pal.list[["H3K4me3"]]) <- c("T=0", "T=3", "T=8")
names(tp.pal.list[["H3K9me3"]]) <- c("T=0", "T=3", "T=8")
names(tp.pal.list[["H3K27ac"]]) <- c("T=0", "T=1", "T=3", "T=8")

st.pal.list <- list()
st.pal.list[["H3K27me3"]] <- c("Leading" = "#000000",
                               "Lagging" = brewer.pal(n=9,name="Purples")[8])
st.pal.list[["H3K4me3"]] <- c("Leading" = "#000000",
                              "Lagging" = brewer.pal(n=9,name="Greens")[8])
st.pal.list[["H3K9me3"]] <- c("Leading" = "#000000",
                              "Lagging" = brewer.pal(n=9,name="Oranges")[8])

st.pal.list[["H3K27ac"]] <- c("Leading" = "#000000",
                              "Lagging" = brewer.pal(n=9,name="Blues")[8])

cond.colors <- c("MCM2-2A" = "#D55E00", "WT" = "#0072B2")

tp.labs <- c("T0" = "T=0", "T1" = "T=1", 
             "T3" = "T=3", "T8" = "T=8")

replicate.labs <- c("r1" = "replicate #1", 
                    "r3" = "replicate #1",
                    "r4" = "replicate #1", 
                    "r2" = "replicate #2",
                    "r5" = "replicate #2")

clone.labs <- c("410" = "WT", "439" = "MCM2-2A #1", 
                "442" = "MCM2-2A #2", "438" = "MCM2-2A #3",
                "423" = "MCM2-Y90A #1",
                "586" = "MCM2-2A-R #1","588" = "MCM2-2A-R #2", 
                "551" = "POLE4-KO #1",
                "552" = "POLE4-KO #2")

tp.rep.labs <- c("rep1\nT=0","rep1\nT=0","rep1\nT=0",
                 "rep2\nT=0","rep2\nT=0",
                 "rep1\nT=1","rep1\nT=1","rep1\nT=1",
                 "rep2\nT=1","rep2\nT=1",
                 "rep1\nT=3","rep1\nT=3","rep1\nT=3",
                 "rep2\nT=3","rep2\nT=3",
                 "rep1\nT=8","rep1\nT=8","rep1\nT=8",
                 "rep2\nT=8","rep2\nT=8")

names(tp.rep.labs) <- c("T0_r1", "T0_r3","T0_r4",
                        "T0_r2", "T0_r5",
                        "T1_r1", "T1_r3","T1_r4",
                        "T1_r2", "T1_r5",
                        "T3_r1", "T3_r3","T3_r4",
                        "T3_r2", "T3_r5",
                        "T8_r1", "T8_r3","T8_r4",
                        "T8_r2", "T8_r5")

clone.conds <- c("410" = "WT", "439" = "MCM2-2A",
                 "442" = "MCM2-2A", "438" = "MCM2-2A",
                 "423" = "MCM2-Y90A",
                 "588" = "MCM2-2A-R",
                 "551" = "POLE4-KO", "552" = "POLE4-KO")

reps <- c("r1" = "r1", "r2" = "r2", "r3" = "r1",
          "r4" = "r1", "r5" = "r2")

tmp.df <- crossing(names(clone.conds), names(reps)) %>%
  as.data.frame()
colnames(tmp.df) <- c("clone", "replicate")
tmp.df$samps <- paste0(tmp.df$clone,"_",tmp.df$replicate)
samp.order <- tmp.df$samps

samp.label <- paste0(clone.labs[tmp.df$clone],"\n", replicate.labs[tmp.df$replicate])
names(samp.label) <- tmp.df$samps

plots.dir <- "plots"
####
load(file.path(chip.data.dir, "ChIP_processed_peaks.RData"))

scar.peaks.list <- list()
scar.peaks.list[["ChIP_H3K27me3"]] <- chip.peaks.list$H3K27me3_410
scar.peaks.list[["ChIP_H3K9me3"]] <- chip.peaks.list$H3K9me3_410
scar.peaks.list[["ChIP_H3K4me3"]] <- chip.peaks.list$H3K4me3_410
scar.peaks.list[["ChIP_H3K27ac"]] <- chip.peaks.list$H3K9me3_410

scar.peaks.list[["ChIP_SUZ12"]] <- chip.peaks.list$SUZ12_410
scar.peaks.list[["ChIP_SUZ12_strict"]] <- chip.peaks.list$SUZ12_410[chip.peaks.list$SUZ12_410$qval >= -log10(0.05)]
scar.peaks.list[["ChIP_H3K27me1"]] <- chip.external.peaks.list$H3K27me1
scar.peaks.list[["ChIP_H3K27me2"]] <- chip.external.peaks.list$H3K27me2
scar.peaks.list[["ChIP_H3K9me2"]] <- chip.external.peaks.list$H3K9me2
scar.peaks.list[["ChIP_H3K36me2"]] <- chip.external.peaks.list$H3K36me2
scar.peaks.list[["ChIP_H3K36me3"]] <- chip.external.peaks.list$H3K36me3
scar.peaks.list[["ChIP_H3K4me1"]] <- chip.external.peaks.list$H3K4me1

scar.peaks.list[["ChIP_H2AK119ub1"]] <- chip.external.peaks.list$H2AK119ub1
scar.peaks.list[["ChIP_H2AK119ub1_strict"]] <- 
  chip.external.peaks.list$H2AK119ub1[chip.external.peaks.list$H2AK119ub1$pval >= -log10(0.05)]

#### Plots of H3K27me3 for MCM2-2A and POLE4-KO

mk <- "H3K27me3"

outfile <- paste0("SCAR_rfd_MCM2", "_",mk,".RData")
load(file.path(scar.data.dir, outfile))

RFD.tmp.gr <- subset(RFD.gr, ((clone  %in% c("410", "439")) & 
                                (timepoint == "T0")) | sample == "OK-seq")

outfile <- paste0("SCAR_rfd_Pole4", "_",mk,".RData")
load(file.path(scar.data.dir, outfile))

RFD.tmp2.gr <- subset(RFD.gr, clone  == "551" & sample != "OK-seq")


RFD.gr <- c(RFD.tmp.gr, RFD.tmp2.gr)
rm(list = c("RFD.tmp.gr", "RFD.tmp2.gr"))
gc()

rfd.dist.chip <- as.data.frame(distanceToNearest(resize(RFD.gr, width = 2, fix = "center"), 
                                                 chip.peaks.list$H3K27me3_410))

elementMetadata(RFD.gr)[["ChIP"]] <- FALSE
elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance <= 100000]])[["ChIP"]] <- TRUE


OK.ext.gr <- subset(RFD.gr, seqnames %in% chroms & IZ & sample == "OK-seq" & ChIP)

OK.ext.gr$break_start <- start(OK.ext.gr)
# exclude overlapping initiation zone (within 200000 bp).
# subset to data within search-space (ok.ext.gr)
OK.ext.gr <- resize(OK.ext.gr,200000,fix="center")


OK.dist <- distanceToNearest(OK.ext.gr)
OK.ext.gr <- OK.ext.gr[-queryHits(subset(OK.dist,
                                         OK.dist@elementMetadata$distance==0))]

rt.labels <- paste0(names(table(OK.ext.gr$RT)), " (n=",
                    as.numeric(table(OK.ext.gr$RT)), ")")

names(rt.labels) <- names(table(OK.ext.gr$RT))

overlap.pairs <- findOverlaps(OK.ext.gr, RFD.gr)
RFD.break.all.gr <- RFD.gr[subjectHits(overlap.pairs)]
RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
RFD.break.all.gr$RT <- OK.ext.gr$RT[queryHits(overlap.pairs)]
RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]


## ===================  Plot a partition plot  ==================================
### Here we define the CPM cutoff
### Higher values usually give higher mean partitions scores
### however also gives more noisy plots since we remove more windows
### play around with it to see what you prefer

CPM.cutoff <- 0.3


### Here just sumarise the mean partion score
### for each sample over the IZs 

RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq", 
                clone %in% c("410","439","551"),
                timepoint == "T0") %>%
  dplyr::group_by(dist,clone,mark,timepoint) %>%  # rank, enh_active
  dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% 
  mutate(condition = clone.conds[clone]) %>%
  as.data.frame()

OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
  dplyr::group_by(dist) %>%  # rank, enh_active
  dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% as.data.frame()

OK.mean.df$clone <- "263"
OK.mean.df$mark <- "OK-seq"
OK.mean.df$timepoint <- "T0"
OK.mean.df$condition <- "OK-seq"


plot.mean.df <- rbind(OK.mean.df, RFD.mean.df)

plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","POLE4-KO", 
                                                                    "OK-seq"))

n.izs <- length(OK.ext.gr)
scar.fig3F.plt <- ggplot(plot.mean.df) + 
  geom_line(stat = "smooth",size = 1.3, aes(x = dist / 1000, y = RFD,
                                            colour=condition),
            method = "gam",se=F, inherit.aes = TRUE) + 
  xlab(paste0("Distance (kb) from initiation zone center")) +
  ylab("Partition or RFD") +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
  scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2), 
                     limits = c(-0.21,0.21)) +
  scale_colour_manual(values=part.cond.cols[levels(plot.mean.df$condition)]) +  guides(colour = guide_legend(title = paste0(mk, "\nSCAR-seq"))) +
  theme_classic(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        aspect.ratio = 0.8) +
  labs(caption = TeX(paste0("$N =",n.izs, "$")),
       family = "Helvetica", size = 5)

ggsave(filename = file.path(plots.dir, paste0("Fig3F_","SCAR","_",mk, "_partition.pdf")),
       plot = scar.fig3F.plt,
       width = 12, height = 10, device = cairo_pdf)


main.part.fig.names <- c("H3K27me3" = "Fig1C",
                            "H3K4me3" = "Fig1D",
                            "H3K27ac" = "Fig1E")

supp.part.rep.fig.names <- c("H3K27me3" = "FigS2A",
                                "H3K4me3" = "FigS3A",
                                "H3K27ac" = "FigS4A")

supp.part.corr.fig.names <- c("H3K27me3" = "FigS2C",
                                "H3K4me3" = "FigS3C",
                                "H3K27ac" = "FigS4C")

supp.part.rt.fig.names <- c("H3K27me3" = "FigS2E",
                                "H3K27ac" = "FigS4E")

supp.part.box.fig.names <- c("H3K27me3" = "FigS2B",
                                 "H3K4me3" = "FigS3B",
                                 "H3K27ac" = "FigS4B")

#### Analysis for all MCM2-2A SCAR-seq H3K27ac/H3K27me3/H3K4me3
prj <- "MCM2"
all.box.stats.df <- c()
all.box.pvals.df <- c()
for(mk in names(main.part.fig.names)) {
  outfile <- paste0("SCAR_rfd_",prj, "_",mk,".RData")
  load(file.path(scar.data.dir, outfile))
  
  print(mk)
  
  pkn <- paste0("ChIP_",mk)
  tmp.gr <- scar.peaks.list[[pkn]]
  tmp.gr <- tmp.gr[seqnames(tmp.gr) %in% chroms]
  seqlevels(tmp.gr) <- chroms
  tmp.gr <- keepSeqlevels(tmp.gr, chroms)
  rfd.dist.chip <- as.data.frame(distanceToNearest(resize(RFD.gr, width = 2, fix = "center"), 
                                                   tmp.gr))
  
  elementMetadata(RFD.gr)[["ChIP"]] <- FALSE
  elementMetadata(RFD.gr)[["ChIP_ovrl"]] <- FALSE
  elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance <= 100000]])[["ChIP"]] <- TRUE
  elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance == 0]])[["ChIP_ovrl"]] <- TRUE
  
  
  OK.ext.gr <- subset(RFD.gr, IZ & sample == "OK-seq" & ChIP)
  
  OK.ext.gr$break_start <- start(OK.ext.gr)
  OK.ext.gr <- resize(OK.ext.gr,200000,fix="center")
  
  
  OK.dist <- distanceToNearest(OK.ext.gr)
  OK.ext.gr <- OK.ext.gr[-queryHits(subset(OK.dist,
                                           OK.dist@elementMetadata$distance==0))]
  
  rt.labels <- paste0(names(table(OK.ext.gr$RT)), " (n=",
                      as.numeric(table(OK.ext.gr$RT)), ")")
  
  names(rt.labels) <- names(table(OK.ext.gr$RT))
  
  overlap.pairs <- findOverlaps(OK.ext.gr, RFD.gr)
  RFD.break.all.gr <- RFD.gr[subjectHits(overlap.pairs)]
  RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
  RFD.break.all.gr$RT <- OK.ext.gr$RT[queryHits(overlap.pairs)]
  RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]
  
  
  ## ===================  Plot a partition plot  ==================================
  ### Here we define the CPM cutoff
  ### Higher values usually give higher mean partitions scores
  ### however also gives more noisy plots since we remove more windows
  ### play around with it to see what you prefer
  
  CPM.cutoff <- 0.3
  
  
  ### Here just sumarise the mean partion score
  ### for each sample over the IZs 
  
  RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
    dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq", 
                  clone %in% c("410","439")) %>%
    dplyr::group_by(dist,clone,mark,timepoint) %>%  # rank, enh_active
    dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
      RFD.raw = mean(RFD.raw, na.rm = T),
      RFD = mean(RFD,na.rm = T)) %>% 
    mutate(condition = ifelse(clone == "410", "WT", "MCM2-2A")) %>%
    as.data.frame()
  
  OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
    dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
    dplyr::group_by(dist) %>%  # rank, enh_active
    dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
      RFD.raw = mean(RFD.raw, na.rm = T),
      RFD = mean(RFD,na.rm = T)) %>% as.data.frame()
  
  OK.mean.df$clone <- "263"
  OK.mean.df$mark <- "OK-seq"
  OK.mean.df$timepoint <- "T0"
  OK.mean.df$condition <- "OK-seq"
  
  
  plot.mean.df <- rbind(OK.mean.df, RFD.mean.df)
  plot.mean.df$timepoint <- factor(gsub("T","T=",plot.mean.df$timepoint),
                                   levels = c("T=0","T=1","T=3","T=8"))
  plot.mean.df$timepoint <- droplevels(plot.mean.df$timepoint)
  
  plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","OK-seq"))
  plot.mean.df$condition <- droplevels(plot.mean.df$condition)
   
  
  ### Here we plot... depending what kind
  ### of factors and how you want to split
  ### the plots you can do one or several
  
  part.cond.cols.tmp <- part.cond.cols[levels(plot.mean.df$condition)]
  tp.alpha.list.tmp <- tp.alpha.list[levels(plot.mean.df$timepoint)]
  
  n.izs <- length(OK.ext.gr)
  fig1.plt <- ggplot(plot.mean.df) + 
    geom_line(stat = "smooth",size = 1.3, aes(x = dist / 1000, y = RFD,
                                              colour=condition, alpha = timepoint),
              method = "gam",se=F, inherit.aes = TRUE) + 
    xlab(paste0("Distance (kb) from initiation zone center")) +
    ylab("Partition or RFD") +
    geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
    scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2), 
                       limits = c(-0.21,0.21)) +
    scale_colour_manual(values=part.cond.cols.tmp) +
    scale_alpha_manual(values=tp.alpha.list.tmp) +
    guides(colour = guide_legend(order = 1, title = paste0(mk, "\nSCAR-seq")),
           alpha = guide_legend(order = 2, title = "")) +
    theme_classic(base_size = 20, base_family = "Helvetica") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          #legend.title = element_blank(),
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 6 * .pt),
          legend.text = element_text(size = 6 * .pt),
          aspect.ratio = 0.8) +
    labs(caption = TeX(paste0("$N =",n.izs, "$")),
         family = "Helvetica", size = 5)
  
  ggsave(filename = file.path(plots.dir, paste0(main.part.fig.names[mk],"_SCAR_",mk, "partition.pdf")),
         plot = fig1.plt,
         width = 12, height = 10, device = cairo_pdf)
  
  
  
  
  RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
    dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq") %>%
    dplyr::group_by(dist,clone,mark,timepoint,replicate) %>%  # rank, enh_active
    dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
      RFD.raw = mean(RFD.raw, na.rm = T),
      RFD = mean(RFD,na.rm = T)) %>% 
    mutate(condition = ifelse(clone == "410", "WT", "MCM2-2A")) %>%
    as.data.frame()
  
  OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
    dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
    dplyr::group_by(dist) %>%  # rank, enh_active
    dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
      RFD.raw = mean(RFD.raw, na.rm = T),
      RFD = mean(RFD,na.rm = T)) %>% as.data.frame()
  
  OK.mean.df$mark <- "OK-seq"
  OK.mean.df$condition <- "OK-seq"
  OK.mean.df$timepoint <- "T0"
  
  plot.mean.df <- c()
  for (cid in unique(RFD.mean.df$clone)) {
    for (rid in unique(RFD.mean.df$replicate)) {
      OK.mean.df$clone <- cid
      OK.mean.df$replicate <- rid
      plot.mean.df <- rbind(plot.mean.df, OK.mean.df)
    }
  }
  plot.mean.df <- rbind(plot.mean.df, RFD.mean.df)
  plot.mean.df$timepoint <- factor(gsub("T","T=",plot.mean.df$timepoint),
                                   levels = c("T=0","T=1","T=3","T=8"))
  
  plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","OK-seq"))
  plot.mean.df$clone <- factor(plot.mean.df$clone, levels = c("410","439","438","442"))
  plot.mean.df$replicate <- factor(reps[as.character(plot.mean.df$replicate)])
  
  max.rfd <- max(abs(plot.mean.df$RFD), na.rm = TRUE) + 0.01
  
  supp.fig.partition.rep.plt  <- 
    ggplot(plot.mean.df) +
    geom_line(stat = "smooth",size = 1, aes(x = dist / 1000, y = RFD,
                                              colour=condition, alpha = timepoint),
              method = "gam",se=F, inherit.aes = TRUE) + 
    xlab(paste0("Distance (kb) from initiation zone center")) +
    ylab("Partition or RFD") +
    geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
    scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2), 
                       limits = c(-max.rfd,max.rfd)) +
    scale_x_continuous(breaks=c(-50, 0, 50), 
                       limits = c(-100,100)) +
    scale_colour_manual(values=part.cond.cols.tmp) +
    scale_alpha_manual(values=tp.alpha.list.tmp) +
    guides(colour = guide_legend(order = 1, title = paste0(mk, "\nSCAR-seq")),
           alpha = guide_legend(order = 2, title = "")) +
    theme_bw(base_size = 20, base_family = "Helvetica") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 6 * .pt),
          legend.text = element_text(size = 6 * .pt),
          aspect.ratio = 0.8) +
    facet_grid(replicate ~ clone, labeller = labeller(.rows = replicate.labs, 
                                                      .cols = clone.labs)) +
    labs(caption = TeX(paste0("$N =",n.izs, "$")),
         family = "Helvetica", size = 5)
  
  ggsave(filename = file.path(plots.dir, paste0(supp.part.rep.fig.names[mk],
                                                "_SCAR_",mk , "_partition_replicate.pdf")),
         plot = supp.fig.partition.rep.plt,
         width = 14, height = 10, device = cairo_pdf)
  
  
  
  RFD.mean.mut.df <- RFD.break.all.gr %>% as.data.frame() %>% 
    dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq",
                  clone != "410") %>%
    dplyr::group_by(dist,clone,mark,timepoint,RT) %>%  # rank, enh_active
    dplyr::summarise(.groups = "keep",
                     RFD.raw = mean(RFD.raw, na.rm = T),
                     RFD = mean(RFD,na.rm = T)) %>% 
    mutate(condition = ifelse(clone == "410", "WT", "MCM2-2A"),
           sample = clone) %>%
    as.data.frame()
  
  OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
    dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
    dplyr::group_by(dist, RT) %>%  # rank, enh_active
    dplyr::summarise(.groups = "keep",
                     RFD.raw = mean(RFD.raw, na.rm = T),
                     RFD = mean(RFD,na.rm = T)) %>% as.data.frame()
  
  OK.mean.df$mark <- "OK-seq"
  OK.mean.df$condition <- "OK-seq"
  OK.mean.df$timepoint <- "T0"
  
  plot.mean.df <- c()
  for (cid in unique(RFD.mean.mut.df$clone)) {
      OK.mean.df$clone <- cid
      OK.mean.df$sample <- cid
      plot.mean.df <- rbind(plot.mean.df, OK.mean.df)
      
      RFD.mean.wt.df <- RFD.break.all.gr %>% as.data.frame() %>% 
        dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq",
                      clone == "410") %>%
        dplyr::group_by(dist,clone,mark,timepoint,RT) %>% 
        dplyr::summarise(.groups = "keep",
                         RFD.raw = mean(RFD.raw, na.rm = T),
                         RFD = mean(RFD,na.rm = T)) %>% 
        mutate(condition = ifelse(clone == "410", "WT", "MCM2-2A")) %>%
        as.data.frame()
      RFD.mean.wt.df$sample <- cid
      plot.mean.df <- rbind(plot.mean.df, RFD.mean.wt.df)
  }
  plot.mean.df <- rbind(plot.mean.df,  RFD.mean.mut.df)
  plot.mean.df$timepoint <- factor(gsub("T","T=",plot.mean.df$timepoint),
                                   levels = c("T=0","T=1","T=3","T=8"))
  
  plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","OK-seq"))
  plot.mean.df$clone <- factor(plot.mean.df$clone, levels = c("410","439","438","442"))
  plot.mean.df$sample <- factor(plot.mean.df$sample, levels = c("439","438","442"))
  
  max.rfd <- max(c(abs(plot.mean.df$RFD),0.15), na.rm = TRUE) + 0.01
  supp.fig.partition.RT.plt <- 

    ggplot(plot.mean.df) +
    geom_line(stat = "smooth",size = 1, aes(x = dist / 1000, y = RFD,
                                            colour=condition, alpha = timepoint),
              method = "gam",se=F, inherit.aes = TRUE) + 
    xlab(paste0("Distance (kb) from initiation zone center")) +
    ylab("Partition or RFD") +
    geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
    scale_y_continuous(breaks=c(-0.15, 0, 0.15), 
                       limits = c(-max.rfd,max.rfd)) +
    scale_x_continuous(breaks=c(-50, 0, 50), 
                       limits = c(-100,100)) +
    scale_colour_manual(values=part.cond.cols.tmp) +
    scale_alpha_manual(values=tp.alpha.list.tmp) +
    guides(colour = guide_legend(order = 1, title = paste0(mk, "\nSCAR-seq")),
           alpha = guide_legend(order = 2, title = "")) +
    theme_bw(base_size = 20, base_family = "Helvetica") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 6 * .pt),
          legend.text = element_text(size = 6 * .pt)) +
    facet_grid(sample ~ RT, labeller = labeller(.rows = clone.labs, .cols = rt.labels)) #+
  
  if (mk != "H3K4me3") {
    ggsave(filename = file.path(plots.dir, paste0(supp.part.rt.fig.names[mk],
                                                  "Rebuttal_SCAR_",mk, "_RT.pdf")),
           plot = supp.fig.partition.RT.plt,
           width = 16, height = 6, device = cairo_pdf)
  }
  
    
  
  RFD.cor.melt.df$clone <- factor(strsplit2(as.character(RFD.cor.melt.df$SCAR),"_")[,1],
                                  levels = names(clone.labs))
  RFD.spear.df$clone <- factor(strsplit2(as.character(RFD.spear.df$SCAR),"_")[,1],
                                  levels = names(clone.labs))
  RFD.cor.melt.df$timepoint <- factor(strsplit2(as.character(RFD.cor.melt.df$SCAR),"_")[,3],
                                  levels = names(tp.labs))
  

  max.part <- 0.7                
  supp.fig.partition.scatter.plt <- RFD.cor.melt.df %>%
    mutate(clone = factor(clone, levels = names(clone.labs))) %>%
    ggplot(aes(x = RFD, y = Partition)) +
    geom_hex(bins=120) + 
    scale_fill_gradientn(colours = (brewer.pal(n=9,name="Blues")[2:8])) +
    geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
    geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
    scale_y_continuous(breaks=c(-0.6, 0, 0.6), 
                       limits = c(-max.part,max.part)) +
    scale_x_continuous(breaks=c(-0.6, 0, 0.6), 
                       limits =  c(-max.part,max.part)) +
    geom_smooth(se=F,method="lm",size=0.3,colour="red") +
    theme_bw(base_size = 20, base_family = "Helvetica") +
    geom_text(data = RFD.spear.df, 
              aes(x = -0.4, y = 0.5, 
                  label = tex.label),
              parse = TRUE, size = 5) +
    facet_grid(clone~timepoint, labeller = labeller(.rows = clone.labs,
                                                    .cols = tp.labs)) +
    theme(panel.spacing = unit(0.3, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 7 * .pt),
      axis.text = element_text(size = 6 * .pt),
      legend.text = element_text(size = 6 * .pt),
      strip.text.x = element_text(size = 5 * .pt),
      strip.text.y = element_text(size = 5 * .pt),
      strip.background = element_rect(colour="black", fill="#e3e3e3"),
      aspect.ratio = 0.8)
  
  ggsave(filename = file.path(plots.dir, paste0(supp.part.corr.fig.names[mk],
                                                "_SCAR_",mk, "_scatter.pdf")),
         plot =   supp.fig.partition.scatter.plt,
         width = 16, height = 12, device = cairo_pdf)
  
  
  ### Boxplots

  mt.ids <- c("439","438")
  bx.breaks <- c(0,3,8)
  tmp.tp.labs <- tp.labs[c(1,3,4)]
  if (mk == "H3K27ac") {
    mt.ids <- c("439","442")
    bx.breaks <- c(0,1,3,8)
    tmp.tp.labs <- tp.labs[c(1,2,3,4)]
  }
  
  OK.ext.gr <- subset(RFD.gr, IZ & sample == "OK-seq")
  
  OK.ext.gr$break_start <- start(OK.ext.gr)
  OK.ext.gr <- resize(OK.ext.gr,200000,fix="center")

  
  OK.dist <- distanceToNearest(OK.ext.gr)
  OK.ext.gr <- OK.ext.gr[-queryHits(subset(OK.dist,
                                           OK.dist@elementMetadata$distance==0))]

  names(rt.labels) <- names(table(OK.ext.gr$RT))
  
  
  overlap.pairs <- findOverlaps(OK.ext.gr, RFD.gr)
  RFD.break.all.gr <- RFD.gr[subjectHits(overlap.pairs)]
  RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
  RFD.break.all.gr$RT <- OK.ext.gr$RT[queryHits(overlap.pairs)]
  RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]
  
  RFD.break.all.gr$direction <- "Upstream"
  RFD.break.all.gr$direction[RFD.break.all.gr$dist > 0] <- "Downstream"
  
  RFD.break.all.gr$RFD.adj <- ifelse(RFD.break.all.gr$direction == "Downstream",
                                     RFD.break.all.gr$RFD, -RFD.break.all.gr$RFD) 
  
  
  
  RFD.break.ok.df <- RFD.break.all.gr %>% as.data.frame() %>% 
    filter(sample == "OK-seq", ChIP)
  
  rownames(RFD.break.ok.df) <- RFD.break.ok.df$names
  reg.to.breaks <- RFD.break.ok.df$break_ID
  names(reg.to.breaks) <- RFD.break.ok.df$names
  
  RFD.ok_F <- RFD.break.ok.df[RFD.break.ok.df$dist>0,] # downstream
  RFD.ok_R <- RFD.break.ok.df[RFD.break.ok.df$dist<0,] # upstream
  
  # highest RFD downstream of initiation zone
  RFD_max <- RFD.ok_F %>% dplyr::group_by(break_ID) %>%
    dplyr::summarise(names = names[which.max(RFD)],
                     RFDedge = max(RFD),
                     dist = dist[which.max(RFD)]) %>% 
    mutate(stand_rep = "Leading") %>%
    as.data.frame()
  RFD_max$seqnames <- strsplit2(RFD_max$names, ":")[,1]
  RFD_max$start <- as.numeric(strsplit2(strsplit2(RFD_max$names, ":")[,2],"-")[,1])
  RFD_max$end <- as.numeric(strsplit2(strsplit2(RFD_max$names, ":")[,2],"-")[,2])
  RFD_max.gr <- makeGRangesFromDataFrame(RFD_max, keep.extra.columns = TRUE)
  
  # lowest RFD upstream of initiation zone
  RFD_min <- RFD.ok_R %>% dplyr::group_by(break_ID) %>%
    dplyr::summarise(names = names[which.min(RFD)],
                     RFDedge = min(RFD),
                     dist = dist[which.min(RFD)]) %>% 
    mutate(stand_rep = "Lagging") %>%
    as.data.frame()
  
  RFD_min$seqnames <- strsplit2(RFD_min$names, ":")[,1]
  RFD_min$start <- as.numeric(strsplit2(strsplit2(RFD_min$names, ":")[,2],"-")[,1])
  RFD_min$end <- as.numeric(strsplit2(strsplit2(RFD_min$names, ":")[,2],"-")[,2])
  RFD_min.gr <- makeGRangesFromDataFrame(RFD_min, keep.extra.columns = TRUE)
  
  brkids <- unique(union(RFD_max.gr$break_ID, RFD_min.gr$break_ID))
  
  tmp.break.all.gr <- subset(tmp.break.all.gr, break_ID %in% brkids)
  
  nearest.max.dist <- as.data.frame(distanceToNearest(RFD.break.all.gr, RFD_max.gr))
  RFD.break.all.gr$nearest.max <- Inf
  RFD.break.all.gr$nearest.max[nearest.max.dist$queryHits] <- nearest.max.dist$distance
  
  nearest.min.dist <- as.data.frame(distanceToNearest(RFD.break.all.gr, RFD_min.gr))
  RFD.break.all.gr$nearest.min <- Inf
  RFD.break.all.gr$nearest.min[nearest.min.dist$queryHits] <- nearest.min.dist$distance
  
  RFD.break.all.gr$distance.extreme <- ifelse(RFD.break.all.gr$direction == "Downstream",
                                              RFD.break.all.gr$nearest.max,
                                              RFD.break.all.gr$nearest.min)
  
  
  RFD.box.df <- RFD.break.all.gr %>% as.data.frame() %>% 
    dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, 
                  distance.extreme != Inf,
                  sample != "OK-seq",
                  mark == mk,
                  abs(dist) >= 10000, abs(dist) <= ((distance.extreme * 2) - 10000)) %>%
    mutate(condition = clone.conds[clone]) %>%
    as.data.frame()
  
  RFD.box.df$timepoint <- factor(RFD.box.df$timepoint,
                                 levels = c("T0","T1","T3","T8"))
  
  RFD.box.df$condition <- factor(RFD.box.df$condition,
                                 levels = c("WT","MCM2-2A"))
  
  tmp.df <- RFD.box.df %>%
    group_by(names, clone, timepoint, condition) %>% 
    summarise(RFD.adj = mean(RFD.adj, na.rm = T)) %>%
    group_by(clone, timepoint, condition) %>%
    mutate(Q3 = quantile(RFD.adj, 0.75), 
           Q1 = quantile(RFD.adj, 0.25), 
           IQR = IQR(RFD.adj)) %>% 
    filter(RFD.adj < (Q3 + 1.5*IQR),  RFD.adj > (Q1 - 1.5 * IQR)) %>%
    mutate(min.RFD = min(RFD.adj, na.rm = TRUE),
           max.RFD = max(RFD.adj, na.rm = TRUE)) %>% select(clone, timepoint,condition,
                                                            min.RFD, max.RFD) %>%
    unique() %>% as.data.frame()
  
  
  box.stats.df <- RFD.box.df %>%
    group_by(names,timepoint, condition, clone) %>% 
    summarise(RFD.adj = mean(RFD.adj, na.rm = TRUE)) %>%
    group_by(timepoint, 
             condition, clone) %>%
    summarise(RFD.mean = mean(RFD.adj, na.rm = TRUE),
              RFD.sd = sd(RFD.adj, na.rm = TRUE),
              RFD.median = median(RFD.adj, na.rm = TRUE),
              n = n()) %>%
    mutate(label = paste0("n=",n),
           clone = factor(clone, levels = c("410",mt.ids)),
           clone.labs2 = factor(clone.labs[as.character(clone)]),
           timepoint = droplevels(timepoint)) %>% as.data.frame()
  
  
  box.stats.df$tp <- as.numeric(gsub("T","",as.character(box.stats.df$timepoint)))
  all.box.stats.df <- rbind(all.box.stats.df, box.stats.df)
  
  RFD.box.df$tp <- as.numeric(gsub("T","",as.character(RFD.box.df$timepoint)))
  
  tmp2.df <- filter(RFD.box.df) %>% 
    group_by(names, clone, condition, timepoint) %>%
    summarise(RFD.adj = mean(RFD.adj, na.rm = TRUE)) %>%
    mutate(clone = factor(clone,levels = c("410",mt.ids)),
           tp = as.numeric(gsub("T","",as.character(timepoint)))) %>%
    mutate(clone.labs2 = factor(clone.labs[as.character(clone)])) %>%
    as.data.frame()  
  
  
  clone.colors.tmp <- c("darkblue","orange","darkred")
  names(clone.colors.tmp) <- clone.labs[c("410",mt.ids)]
  

  y.lim <- max(abs(c(tmp.df$min.RFD, tmp.df$max.RFD)))
  new.partition.boxplot.plt <- ggplot(tmp2.df) +
    geom_boxplot(aes(x = tp, y = RFD.adj, fill = clone.labs2, 
                     group = interaction(tp, clone.labs2)),
                 outlier.shape = NA, width = 0.8,
                 notch = 1, lwd = 1, alpha = 0.8) + 
    geom_line(data = box.stats.df,aes(x = tp, y = RFD.median, color = clone.labs2),
              size = 0.7, position=position_dodge(1), 
              linetype = "dashed") +
    geom_point(data = box.stats.df, aes(x = tp, y = RFD.median, 
                                        color = clone.labs2),
               size = 1.3, position=position_dodge(1)) +
    scale_fill_manual(values=clone.colors.tmp) +
    scale_color_manual(values=clone.colors.tmp) +
    ylim(-y.lim  - 0.015, y.lim + 0.015) +
    ylab("Partition") + xlab("") +
    theme_bw(base_size = 20, base_family = "Helvetica") + 
    scale_x_continuous(breaks = bx.breaks,
                       labels=tmp.tp.labs) +
    guides(fill = guide_legend(title = mk),
           color = "none") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title = element_text(size = 7 * .pt),
          axis.text = element_text(size = 6 * .pt),
          legend.text = element_text(size = 6 * .pt),
          panel.spacing = unit(0.35, "lines"),
          strip.text.x = element_text(size = 6 * .pt),
          strip.text.y = element_text(size = 6 * .pt),
          strip.background = element_rect(colour="black", fill="#e3e3e3")) 
  
  
  ggsave(filename = file.path(plots.dir, paste0(supp.part.box.fig.names[mk],
                                                "_SCAR_",mk,"_partition_box.pdf")),
         plot =   new.partition.boxplot.plt ,
         width = 14, height = 8, device = cairo_pdf)
  
  
  box.pvals.df <- compare_means(data = tmp2.df,
                      formula = RFD.adj ~ clone.labs2, 
                      group.by = "timepoint",
                      method = "wilcox.test") %>% 
    mutate(mark = mk) %>%
    as.data.frame()
  
  
  box.es.df <- c()
  for(tpi in levels(droplevels(tmp2.df$timepoint))) {
    tmp.es.df <- wilcox_effsize(data = tmp2.df[tmp2.df$timepoint == tpi,],
                                formula = RFD.adj ~ clone.labs2,
                                method = "wilcox.test") %>%
      mutate(timepoint = tpi) %>% as.data.frame()
    box.es.df <- rbind(box.es.df, tmp.es.df)
  }
  
  box.wt.df <- full_join(box.es.df, box.pvals.df)
  
  all.box.pvals.df <- rbind(all.box.pvals.df, box.wt.df)
}


write.xlsx(all.box.pvals.df, 
           file = file.path(scar.data.dir,"SCAR_boxplots_stats.xlsx"))

mark.linetype.list <- c("SUZ12" = "solid", "H3K27me3" = "dashed")

prj <- "Crosslinked"
SCAR.all.gr <- GRanges()
i <- 1
for (mk in c("H3K27me3", "SUZ12")) {
  outfile <- paste0("SCAR_rfd_",prj, "_",mk,".RData")
  load(file.path(scar.data.dir, outfile))
  
  print(mk)
  RFD.gr$num <- i
  for (pkn in paste0("ChIP_",c("H3K27me3","SUZ12_strict"))) {
    tmp.gr <- scar.peaks.list[[pkn]]
    pkn2 <- paste0(pkn,"_ovrl")
    tmp.gr <- tmp.gr[seqnames(tmp.gr) %in% chroms]
    seqlevels(tmp.gr) <- chroms
    tmp.gr <- keepSeqlevels(tmp.gr, chroms)
    rfd.dist.chip <- as.data.frame(distanceToNearest(resize(RFD.gr, width = 2, fix = "center"), 
                                                     tmp.gr))
    
    elementMetadata(RFD.gr)[[pkn]] <- FALSE
    elementMetadata(RFD.gr)[[pkn2]] <- FALSE
    elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance <= 100000]])[[pkn]] <- TRUE
    elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance == 0]])[[pkn2]] <- TRUE
  }
  
  SCAR.all.gr <- c(SCAR.all.gr, RFD.gr)
  i <- i+1
}

SCAR.all.gr <- subset(SCAR.all.gr, !((sample == "OK-seq") & (num == 2)))
OK.ext.gr <- subset(SCAR.all.gr, IZ & sample == "OK-seq" & ChIP_H3K27me3 & ChIP_SUZ12_strict)

OK.ext.gr$break_start <- start(OK.ext.gr)
OK.ext.gr <- resize(OK.ext.gr,200000,fix="center")


OK.dist <- distanceToNearest(OK.ext.gr)
OK.ext.gr <- OK.ext.gr[-queryHits(subset(OK.dist,
                                         OK.dist@elementMetadata$distance==0))]
rt.labels <- paste0(names(table(OK.ext.gr$RT)), " (n=",
                    as.numeric(table(OK.ext.gr$RT)), ")")

names(rt.labels) <- names(table(OK.ext.gr$RT))


overlap.pairs <- findOverlaps(OK.ext.gr, SCAR.all.gr)
RFD.break.all.gr <- SCAR.all.gr[subjectHits(overlap.pairs)]
RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
RFD.break.all.gr$RT <- OK.ext.gr$RT[queryHits(overlap.pairs)]
RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]


## ===================  Plot a partition plot  ==================================
### Here we define the CPM cutoff
### Higher values usually give higher mean partitions scores
### however also gives more noisy plots since we remove more windows
### play around with it to see what you prefer

CPM.cutoff <- 0.3

RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq", 
                clone %in% c("410","439")) %>%
  dplyr::group_by(dist,clone,mark) %>%  # rank, enh_active
  dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% 
  mutate(condition = ifelse(clone == "410", "WT", "MCM2-2A")) %>%
  as.data.frame()

OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
  dplyr::group_by(dist) %>%  # rank, enh_active
  dplyr::summarise(#RFD_sd = sd(RFD,na.rm = T),
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% as.data.frame()

OK.mean.df$clone <- "263"
OK.mean.df$mark <- "SUZ12"
OK.mean.df$condition <- "OK-seq"


plot.mean.df <- rbind(OK.mean.df, RFD.mean.df)
plot.mean.df$mark <- factor(plot.mean.df$mark, levels = c("H3K27me3", "SUZ12", "OK-seq"))
plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","OK-seq"))


mx.lim <- max(abs(plot.mean.df$RFD)) + 0.01

n.izs <- length(OK.ext.gr)
scar.suz12.strict.plt <- ggplot(plot.mean.df) + 
  geom_line(stat = "smooth",size = 1.3, aes(x = dist / 1000, y = RFD,
                                            colour=condition, linetype = mark),
            method = "gam",se=F, inherit.aes = TRUE) + 
  xlab(paste0("Distance (kb) from initiation zone center")) +
  ylab("Partition or RFD") +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
  scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2), 
                     limits = c(-mx.lim ,mx.lim)) +
  scale_colour_manual(values=part.cond.cols[c("OK-seq","WT","MCM2-2A")]) +
  scale_linetype_manual(values = mark.linetype.list) +
  guides(colour = guide_legend(title = ""),
         linetype = guide_legend(order = 2, title = "SCAR-seq",
                                 override.aes = list(linetype = c("H3K27me3" = 6,"SUZ12" = 1)))) +
  theme_classic(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        aspect.ratio = 0.8) +
  labs(caption = TeX(paste0("$N =",n.izs, "$")),
       family = "Helvetica", size = 5)

ggsave(filename = file.path(plots.dir, paste0("Fig2E_","SCAR_SUZ12","partition.pdf")),
       plot = scar.suz12.strict.plt,
       width = 12, height = 10, device = cairo_pdf)



RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq") %>%
  dplyr::group_by(dist,clone,mark,replicate) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% 
  mutate(condition = ifelse(clone == "410", "WT", "MCM2-2A")) %>%
  as.data.frame()

plot.mean.df <- c()

OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
  dplyr::group_by(dist) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% as.data.frame()

for (cid in unique(RFD.mean.df$clone)) {
  for (rid in unique(RFD.mean.df$replicate)) {
    OK.mean.df$clone <- cid
    OK.mean.df$mark <- "SUZ12"
    OK.mean.df$condition <- "OK-seq"
    OK.mean.df$replicate <- rid
    plot.mean.df <- rbind(plot.mean.df, OK.mean.df)
  }
}


plot.mean.df <- rbind(plot.mean.df, RFD.mean.df)

plot.mean.df$mark <- factor(plot.mean.df$mark, levels = c("H3K27me3", "SUZ12", "OK-seq"))
plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","OK-seq"))
plot.mean.df$clone <- factor(plot.mean.df$clone, levels = names(clone.labs))
plot.mean.df$replicate<- factor(plot.mean.df$replicate)


max.rfd <- max(abs(plot.mean.df$RFD), na.rm = TRUE) + 0.01

supp.fig.suz12.partition.rep.plt  <- 
  ggplot(plot.mean.df) +
  geom_line(stat = "smooth",size = 1, aes(x = dist / 1000, y = RFD,
                                          colour = condition, 
                                          linetype = mark),
            method = "gam",se=F, inherit.aes = TRUE) + 
  xlab(paste0("Distance (kb) from initiation zone center")) +
  ylab("Partition or RFD") +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
  scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2)) +
  scale_x_continuous(breaks=c(-50, 0, 50), 
                     limits = c(-100,100)) +
  scale_colour_manual(values=part.cond.cols[c("OK-seq","WT","MCM2-2A")]) +
  scale_linetype_manual(values = mark.linetype.list) +
  guides(colour = guide_legend(title = ""),
         linetype = guide_legend(order = 2, title = "SCAR-seq",
                                 override.aes = list(linetype = c("H3K27me3" = 6,"SUZ12" = 1)))) +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        aspect.ratio = 0.8) +
  facet_grid(replicate ~ clone, labeller = labeller(.rows = replicate.labs, 
                                                    .cols = clone.labs)) +
  labs(caption = TeX(paste0("$N =",n.izs, "$")),
       family = "Helvetica", size = 5)

ggsave(filename = file.path(plots.dir, paste0("FigS8D_SCAR_SUZ12_partition_replicate.pdf")),
       plot = supp.fig.suz12.partition.rep.plt,
       width = 14, height = 10, device = cairo_pdf)

### Partition analyses for Rescue samples of H3K27me3/H4K20me0

prj <- "Rescue"
SCAR.all.gr <- GRanges()
i <- 1
for (mk in c("H3K27me3", "H4K20me0")) {
  outfile <- paste0("SCAR_rfd_",prj, "_",mk,".RData")
  load(file.path("data/SCAR_seq", outfile))
  
  print(mk)
  RFD.gr$num <- i
  for (pkn in paste0("ChIP_",c("H3K27me3"))) {
    tmp.gr <- scar.peaks.list[[pkn]]
    pkn2 <- paste0(pkn,"_ovrl")
    tmp.gr <- tmp.gr[seqnames(tmp.gr) %in% chroms]
    seqlevels(tmp.gr) <- chroms
    tmp.gr <- keepSeqlevels(tmp.gr, chroms)
    rfd.dist.chip <- as.data.frame(distanceToNearest(resize(RFD.gr, width = 2, fix = "center"), 
                                                     tmp.gr))
    
    elementMetadata(RFD.gr)[[pkn]] <- FALSE
    elementMetadata(RFD.gr)[[pkn2]] <- FALSE
    elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance <= 100000]])[[pkn]] <- TRUE
    elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance == 0]])[[pkn2]] <- TRUE
  }
  
  SCAR.all.gr <- c(SCAR.all.gr, RFD.gr)
  i <- i+1
}

SCAR.all.gr <- subset(SCAR.all.gr, !((sample == "OK-seq") & (num == 2)))
OK.ext.gr <- subset(SCAR.all.gr, IZ & sample == "OK-seq" & ChIP_H3K27me3)

OK.ext.gr$break_start <- start(OK.ext.gr)
OK.ext.gr <- resize(OK.ext.gr,200000,fix="center")


OK.dist <- distanceToNearest(OK.ext.gr)
OK.ext.gr <- OK.ext.gr[-queryHits(subset(OK.dist,
                                         OK.dist@elementMetadata$distance==0))]
rt.labels <- paste0(names(table(OK.ext.gr$RT)), " (n=",
                    as.numeric(table(OK.ext.gr$RT)), ")")

names(rt.labels) <- names(table(OK.ext.gr$RT))


overlap.pairs <- findOverlaps(OK.ext.gr, SCAR.all.gr)
RFD.break.all.gr <- SCAR.all.gr[subjectHits(overlap.pairs)]
RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
RFD.break.all.gr$RT <- OK.ext.gr$RT[queryHits(overlap.pairs)]
RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]


## ===================  Plot a partition plot  ==================================
### Here we define the CPM cutoff
### Higher values usually give higher mean partitions scores
### however also gives more noisy plots since we remove more windows
### play around with it to see what you prefer


RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq",
                clone %in% c("410","442","588")) %>%
  dplyr::group_by(dist,clone,mark,replicate) %>%
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% 
  mutate(condition = clone.conds[clone]) %>%
  as.data.frame()

plot.mean.df <- c()

OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
  dplyr::group_by(dist) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% as.data.frame()

for (cid in unique(RFD.mean.df$clone)) {
  for (rid in unique(RFD.mean.df$replicate)) {
    OK.mean.df$clone <- cid
    OK.mean.df$mark <- "OK-seq"
    OK.mean.df$condition <- "OK-seq"
    OK.mean.df$replicate <- rid
    plot.mean.df <- rbind(plot.mean.df, OK.mean.df)
  }
}


plot.mean.df <- rbind(plot.mean.df, RFD.mean.df)

plot.mean.df$mark <- factor(plot.mean.df$mark, levels = c("H3K27me3", "H4K20me0", "OK-seq"))
plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","MCM2-2A-R","OK-seq"))
plot.mean.df$clone <- factor(plot.mean.df$clone, levels = c("410","442","588"))
plot.mean.df$replicate <- factor(plot.mean.df$replicate)


mark.linetype.list <- c("H3K27me3" = "dashed","H4K20me0" = "dotted",
                        "OK-seq" = "solid")

n.izs <- length(OK.ext.gr)

supp.fig.rescue.partition.rep.plt  <- 
  ggplot(plot.mean.df) +
  geom_line(stat = "smooth",size = 1, aes(x = dist / 1000, y = RFD,
                                          colour = condition, 
                                          linetype = mark),
            method = "gam",se=F, inherit.aes = TRUE) + 
  xlab(paste0("Distance (kb) from initiation zone center")) +
  ylab("Partition or RFD") +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
  scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2)) +
  scale_x_continuous(breaks=c(-50, 0, 50), 
                     limits = c(-100,100)) +
  scale_colour_manual(values=part.cond.cols[c("OK-seq","WT","MCM2-2A", "MCM2-2A-R")]) +
  scale_linetype_manual(values = mark.linetype.list) +
  guides(colour = guide_legend(title = ""),
         linetype = guide_legend(order = 2, title = "SCAR-seq")) +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        aspect.ratio = 0.8) +
  facet_grid(replicate ~ clone, labeller = labeller(.rows = replicate.labs, 
                                                    .cols = clone.labs)) +
  labs(caption = TeX(paste0("$N =",n.izs, "$")),
       family = "Helvetica", size = 5)

ggsave(filename = file.path(plots.dir, paste0("FigS13A_SCAR_Rescue_partition_replicate.pdf")),
       plot = supp.fig.rescue.partition.rep.plt,
       width = 14, height = 10, device = cairo_pdf)



RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq",
                clone %in% c("410","442","423")) %>%
  dplyr::group_by(dist,clone,mark,replicate) %>% 
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% 
  mutate(condition = clone.conds[clone]) %>%
  as.data.frame()

plot.mean.df <- c()

OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
  dplyr::group_by(dist) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% as.data.frame()

for (cid in unique(RFD.mean.df$clone)) {
  for (rid in unique(RFD.mean.df$replicate)) {
    OK.mean.df$clone <- cid
    OK.mean.df$mark <- "OK-seq"
    OK.mean.df$condition <- "OK-seq"
    OK.mean.df$replicate <- rid
    plot.mean.df <- rbind(plot.mean.df, OK.mean.df)
  }
}


plot.mean.df <- rbind(plot.mean.df, RFD.mean.df)

plot.mean.df$mark <- factor(plot.mean.df$mark, levels = c("H3K27me3", "H4K20me0", "OK-seq"))
plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","MCM2-Y90A","OK-seq"))
plot.mean.df$clone <- factor(plot.mean.df$clone, levels = c("410","442","423"))
plot.mean.df$replicate <- factor(plot.mean.df$replicate)


mark.linetype.list <- c("H3K27me3" = "dashed","H4K20me0" = "dotted",
                        "OK-seq" = "solid")

n.izs <- length(OK.ext.gr)


rescue.plot.mean.df <- plot.mean.df

part.cond.cols["MCM2-2A"] <- brewer.pal(8,"Oranges")[6]
part.cond.cols["MCM2-Y90A"] <- brewer.pal(8,"Oranges")[3]


### Single mutant partition 
y90a.fig.rescue.partition.rep.plt  <- 
  filter(rescue.plot.mean.df,clone %in% c("410","442","423"),
         mark %in% c("OK-seq","H3K27me3")) %>%
  group_by(dist, clone, mark, condition) %>%
  summarise(RFD = mean(RFD)) %>%
  mutate(condition = factor(condition,
                            levels = c("WT","MCM2-2A","MCM2-Y90A","OK-seq"))) %>%
  ggplot() +
  geom_line(stat = "smooth",size = 1, aes(x = dist / 1000, y = RFD,
                                          colour = condition, 
                                          linetype = mark),
            method = "gam",se=F, inherit.aes = TRUE) + 
  xlab(paste0("Distance (kb) from initiation zone center")) +
  ylab("Partition or RFD") +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
  scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2)) +
  scale_x_continuous(breaks=c(-50, 0, 50), 
                     limits = c(-100,100)) +
  scale_colour_manual(values=part.cond.cols[c("WT","MCM2-2A","MCM2-Y90A","OK-seq")]) +
  scale_linetype_manual(values = tmp.mark.linetype.list[c("OK-seq","H3K27me3")]) +
  guides(colour = guide_legend(title = "SCAR-seq\nH3K27me3"),
         linetype = "none") +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        aspect.ratio = 0.8) +
  labs(caption = TeX(paste0("$N =",n.izs, "$")),
       family = "Helvetica", size = 5)

ggsave(filename = file.path(plots.dir, paste0("FigR2D_SCAR_y90A_partition_mean_merged.pdf")),
       plot = y90a.fig.rescue.partition.rep.plt,
       width = 14, height = 12, device = cairo_pdf)


### Partition analyses for POLE4-KO samples of H3K27me3/H4K20me0

prj <- "Pole4"
SCAR.all.gr <- GRanges()
i <- 1
for (mk in c("H3K27me3", "H4K20me0")) {
  outfile <- paste0("SCAR_rfd_",prj, "_",mk,".RData")
  load(file.path(scar.data.dir, outfile))
  
  print(mk)
  RFD.gr$num <- i
  for (pkn in paste0("ChIP_",c("H3K27me3"))) {
    tmp.gr <- scar.peaks.list[[pkn]]
    pkn2 <- paste0(pkn,"_ovrl")
    tmp.gr <- tmp.gr[seqnames(tmp.gr) %in% chroms]
    seqlevels(tmp.gr) <- chroms
    tmp.gr <- keepSeqlevels(tmp.gr, chroms)
    rfd.dist.chip <- as.data.frame(distanceToNearest(resize(RFD.gr, width = 2, fix = "center"), 
                                                     tmp.gr))
    
    elementMetadata(RFD.gr)[[pkn]] <- FALSE
    elementMetadata(RFD.gr)[[pkn2]] <- FALSE
    elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance <= 100000]])[[pkn]] <- TRUE
    elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance == 0]])[[pkn2]] <- TRUE
  }
  
  SCAR.all.gr <- c(SCAR.all.gr, RFD.gr)
  i <- i+1
}

SCAR.all.gr <- subset(SCAR.all.gr, !((sample == "OK-seq") & (num == 2)))
OK.ext.gr <- subset(SCAR.all.gr, IZ & sample == "OK-seq" & ChIP_H3K27me3)

OK.ext.gr$break_start <- start(OK.ext.gr)
OK.ext.gr <- resize(OK.ext.gr,200000,fix="center")

OK.dist <- distanceToNearest(OK.ext.gr)
OK.ext.gr <- OK.ext.gr[-queryHits(subset(OK.dist,
                                         OK.dist@elementMetadata$distance==0))]
rt.labels <- paste0(names(table(OK.ext.gr$RT)), " (n=",
                    as.numeric(table(OK.ext.gr$RT)), ")")

names(rt.labels) <- names(table(OK.ext.gr$RT))


overlap.pairs <- findOverlaps(OK.ext.gr, SCAR.all.gr)
RFD.break.all.gr <- SCAR.all.gr[subjectHits(overlap.pairs)]
RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
RFD.break.all.gr$RT <- OK.ext.gr$RT[queryHits(overlap.pairs)]
RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]


RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq",
                clone %in% c("410","551","552")) %>%
  dplyr::group_by(dist,clone,mark,replicate) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% 
  mutate(condition = clone.conds[clone]) %>%
  as.data.frame()

plot.mean.df <- c()

OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
  dplyr::group_by(dist) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% as.data.frame()

for (cid in unique(RFD.mean.df$clone)) {
  for (rid in unique(RFD.mean.df$replicate)) {
    OK.mean.df$clone <- cid
    OK.mean.df$mark <- "OK-seq"
    OK.mean.df$condition <- "OK-seq"
    OK.mean.df$replicate <- rid
    plot.mean.df <- rbind(plot.mean.df, OK.mean.df)
  }
}


plot.mean.df <- rbind(plot.mean.df, RFD.mean.df)

plot.mean.df$mark <- factor(plot.mean.df$mark, levels = c("H3K27me3", "H4K20me0", "OK-seq"))
plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","POLE4-KO","OK-seq"))
plot.mean.df$clone <- factor(plot.mean.df$clone, levels = c("410","551","552"))
plot.mean.df$replicate <- factor(plot.mean.df$replicate)


mark.linetype.list <- c("H3K27me3" = "dashed","H4K20me0" = "dotted",
                        "OK-seq" = "solid")

n.izs <- length(OK.ext.gr)

supp.fig.pole4.partition.rep.plt  <- 
  ggplot(plot.mean.df) +
  geom_line(stat = "smooth",size = 1, aes(x = dist / 1000, y = RFD,
                                          colour = condition, 
                                          linetype = mark),
            method = "gam",se=F, inherit.aes = TRUE) + 
  xlab(paste0("Distance (kb) from initiation zone center")) +
  ylab("Partition or RFD") +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
  scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2)) +
  scale_x_continuous(breaks=c(-50, 0, 50), 
                     limits = c(-100,100)) +
  scale_colour_manual(values=part.cond.cols[c("OK-seq","WT","POLE4-KO")]) +
  scale_linetype_manual(values = mark.linetype.list) +
  guides(colour = guide_legend(title = ""),
         linetype = guide_legend(order = 2, title = "SCAR-seq",
                                 override.aes = list(linetype = c("H3K27me3" = 6,"H4K20me0" = 3, "OK-seq" = 1)))) +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt),
        aspect.ratio = 0.8) +
  facet_grid(replicate ~ clone, labeller = labeller(.rows = replicate.labs, 
                                                    .cols = clone.labs)) +
  labs(caption = TeX(paste0("$N =",n.izs, "$")),
       family = "Helvetica", size = 5)

ggsave(filename = file.path(plots.dir, paste0("FigS12C_SCAR_Pole4_partition_replicate.pdf")),
       plot = supp.fig.pole4.partition.rep.plt,
       width = 14, height = 10, device = cairo_pdf)


##=========== Plots for stranded inputs ==============================##

prj <- "Input"
outfile <- paste0("SCAR_rfd_",prj, "_","input",".RData")
load(file.path(scar.data.dir, outfile))

OK.ext.gr <- subset(RFD.gr, IZ & sample == "OK-seq")

OK.ext.gr$break_start <- start(OK.ext.gr)
OK.ext.gr <- resize(OK.ext.gr,200000,fix="center")

OK.dist <- distanceToNearest(OK.ext.gr)
OK.ext.gr <- OK.ext.gr[-queryHits(subset(OK.dist,
                                         OK.dist@elementMetadata$distance==0))]
rt.labels <- paste0(names(table(OK.ext.gr$RT)), " (n=",
                    as.numeric(table(OK.ext.gr$RT)), ")")

names(rt.labels) <- names(table(OK.ext.gr$RT))

CPM.cutoff <- 0.03
RFD.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq",
                clone %in% c("439","442"), mark == "input") %>%
  dplyr::group_by(dist,clone,timepoint,mark) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% 
  mutate(condition = clone.conds[clone],
         sample = clone) %>%
  as.data.frame()

OK.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample == "OK-seq") %>%
  dplyr::group_by(dist) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% as.data.frame()

WT.mean.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq",
                clone %in% c("410"), mark == "input") %>%
  dplyr::group_by(dist,clone,timepoint,mark) %>%  
  dplyr::summarise(
    RFD.raw = mean(RFD.raw, na.rm = T),
    RFD = mean(RFD,na.rm = T)) %>% 
  mutate(condition = clone.conds[clone]) %>%
  as.data.frame()

plot.mean.df <- c()

for (sid in unique(RFD.mean.df$sample)) {
  WT.mean.df$sample <- sid
  plot.mean.df <- rbind(plot.mean.df,  WT.mean.df)
}

for (tid in unique(RFD.mean.df$timepoint)) {
  for (cid in unique(RFD.mean.df$clone)) {
    OK.mean.df$clone <- cid
    OK.mean.df$mark <- "OK-seq"
    OK.mean.df$condition <- "OK-seq"
    OK.mean.df$timepoint <- tid
    OK.mean.df$sample <- cid
    plot.mean.df <- rbind(plot.mean.df, OK.mean.df) 
  }
}



plot.mean.df <- rbind(plot.mean.df, RFD.mean.df)

plot.mean.df$condition <- factor(plot.mean.df$condition, levels = c("WT","MCM2-2A","OK-seq"))
plot.mean.df$clone <- factor(plot.mean.df$clone, levels = c("410","439","442"))
plot.mean.df$sample <- factor(plot.mean.df$sample, levels = c("439","442"))
plot.mean.df$timepoint <- factor(plot.mean.df$timepoint)


n.izs <- length(OK.ext.gr)
supp.fig.inputs.partition.rep.plt <- ggplot(plot.mean.df) +
  geom_line(stat = "smooth",size = 1, aes(x = dist / 1000, y = RFD,
                                          colour = condition),
            method = "gam",se=F, inherit.aes = TRUE) + 
  xlab(paste0("Distance (kb) from initiation zone center")) +
  ylab("Partition or RFD") +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
  scale_y_continuous(breaks=c(-0.2,-0.1, 0, 0.1,0.2)) +
  scale_x_continuous(breaks=c(-50, 0, 50), 
                     limits = c(-100,100)) +
  scale_colour_manual(values=part.cond.cols[c("OK-seq","WT","MCM2-2A")]) +
  scale_linetype_manual(values = mark.linetype.list) +
  guides(colour = guide_legend(title = "SCAR-seq\nStranded Inputs")) +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 6 * .pt),
        legend.text = element_text(size = 6 * .pt)) +
  facet_grid(sample ~ timepoint, labeller = labeller(.rows = clone.labs, 
                                                     .cols = tp.labs)) +
  labs(caption = TeX(paste0("$N =",n.izs, "$")),
       family = "Helvetica", size = 5)


ggsave(filename = file.path(plots.dir, paste0("SCAR_FigS4G_Inputs_partition_replicate.pdf")),
       plot = supp.fig.inputs.partition.rep.plt,
       width = 16, height = 10, device = cairo_pdf)




overlap.pairs <- findOverlaps(OK.ext.gr, RFD.gr)
RFD.break.all.gr <- RFD.gr[subjectHits(overlap.pairs)]
RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
RFD.break.all.gr$RT <- OK.ext.gr$RT[queryHits(overlap.pairs)]
RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]

RFD.break.all.gr$direction <- "Upstream"
RFD.break.all.gr$direction[RFD.break.all.gr$dist > 0] <- "Downstream"

RFD.break.all.gr$RFD.adj <- ifelse(RFD.break.all.gr$direction == "Downstream",
                                   RFD.break.all.gr$RFD, -RFD.break.all.gr$RFD) 



CPM.cutoff <- 1
RFD.box.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  dplyr::filter((F.cpm + R.cpm) >= CPM.cutoff, sample != "OK-seq",
                clone %in% c("410","439","442"), mark == "input",
                abs(dist) >= 15000, abs(dist) <= 45000) %>%
  mutate(condition = clone.conds[clone],
         sample = clone) %>%
  as.data.frame()

RFD.box.df$timepoint <- factor(RFD.box.df$timepoint,
                               levels = c("T0","T1","T3","T8"))

RFD.box.df$condition <- factor(RFD.box.df$condition,
                               levels = c("WT","MCM2-2A"))


tmp.df <- RFD.box.df %>%
  group_by(clone, timepoint, condition) %>% mutate(Q3 = quantile(RFD.adj, 0.75), 
                                         Q1 = quantile(RFD.adj, 0.25), 
                                         IQR = IQR(RFD.adj)) %>%
  filter(RFD < (Q3 + 1.5*IQR),  RFD > (Q1 - 1.5 * IQR)) %>%
  mutate(min.RFD = min(RFD.adj, na.rm = TRUE),
         max.RFD = max(RFD.adj, na.rm = TRUE)) %>% select(clone, timepoint,condition, 
                                                      min.RFD, max.RFD) %>%
  unique() %>% as.data.frame()


box.stats.df <- RFD.box.df %>%
  group_by(names,timepoint, condition, clone) %>% 
  summarise(RFD.adj = mean(RFD.adj, na.rm = TRUE)) %>% group_by(timepoint, condition, clone) %>%
  tally() %>%
  mutate(label = paste0("n=",n)) %>% as.data.frame()


inputs.box.439.plt <- filter(RFD.box.df, clone != "442") %>% group_by(names, clone, condition, timepoint) %>%
  summarise(RFD.adj = mean(RFD.adj, na.rm = TRUE)) %>%
  ggplot(aes(x = timepoint, y = RFD.adj, fill = condition)) +
  geom_boxplot(aes(fill=condition),
               outlier.shape = NA,
               notch = 1, lwd = 0.8) + 
  ggtitle("MCM2-2A #1") +
  scale_fill_manual(values=part.cond.cols[c("WT","MCM2-2A")]) +
  geom_text(data = filter(box.stats.df, clone != "442"), 
            aes(x = timepoint, y = min(tmp.df$min.RFD) - 0.05,
                label = label, group = condition),
            position=position_dodge(width=1)) +
  ylab("Partition") + xlab("") +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  scale_x_discrete(labels=tp.labs) +
  labs(fill = "", colour = "Strand", alpha = "Strand") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 5 * .pt),
        legend.text = element_text(size = 6 * .pt),
        panel.spacing = unit(0.35, "lines"),
        strip.text.x = element_text(size = 5 * .pt),
        strip.text.y = element_text(size = 5 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3")) +
  stat_compare_means(aes(group = condition, label = as.character(scientific_10x(..p.adj..))), 
                     method = "wilcox.test", size = 4.5, parse = TRUE, 
                     family = "Helvetica", fontface = "bold", label.y = max(tmp.df$max.RFD) + 0.05,
                     paired = FALSE, na.rm = TRUE)


inputs.box.442.plt <- filter(RFD.box.df, clone != "439") %>% group_by(names, clone, condition, timepoint) %>%
  summarise(RFD.adj = mean(RFD.adj, na.rm = TRUE)) %>%
ggplot(aes(x = timepoint, y = RFD.adj, fill = condition)) +
  geom_boxplot(aes(fill=condition),
               outlier.shape = NA,
               notch = 1, lwd = 0.8) + 
  ggtitle("MCM2-2A #2") +
  scale_fill_manual(values=part.cond.cols[c("WT","MCM2-2A")]) +
  geom_text(data = filter(box.stats.df, clone != "439"), 
            aes(x = timepoint, y = min(tmp.df$min.RFD) - 0.05,
                                     label = label, group = condition),
            position=position_dodge(width=1)) +
  ylab("Partition") + xlab("") +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  scale_x_discrete(labels=tp.labs) +
  labs(fill = "", colour = "Strand", alpha = "Strand") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 5 * .pt),
        legend.text = element_text(size = 6 * .pt),
        panel.spacing = unit(0.35, "lines"),
        strip.text.x = element_text(size = 5 * .pt),
        strip.text.y = element_text(size = 5 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3")) +
  stat_compare_means(aes(group = condition, label = as.character(scientific_10x(..p.adj..))), 
                     method = "wilcox.test", size = 4.5, parse = TRUE, 
                     family = "Helvetica", fontface = "bold", label.y = max(tmp.df$max.RFD) + 0.05,
                     paired = FALSE, na.rm = TRUE)

box.plt <- (inputs.box.439.plt | inputs.box.442.plt) + 
  plot_layout(guides = "collect")
ggsave(filename = file.path(plots.dir,paste0("FigS4H_SCAR_","Inputs", "_boxplots_cond_CPM",CPM.cutoff,
                         ".pdf")),
       plot =   box.plt,
       width = 18, height = 10, device = cairo_pdf)




RFD.cor.melt.df$clone <- factor(strsplit2(as.character(RFD.cor.melt.df$SCAR),"_")[,1],
                                levels = names(clone.labs))
RFD.spear.df$clone <- factor(strsplit2(as.character(RFD.spear.df$SCAR),"_")[,1],
                             levels = names(clone.labs))
RFD.cor.melt.df$timepoint <- factor(strsplit2(as.character(RFD.cor.melt.df$SCAR),"_")[,3],
                                    levels = names(tp.labs))
RFD.spear.df$timepoint <- factor(strsplit2(as.character(RFD.spear.df$SCAR),"_")[,3],
                                 levels = names(tp.labs))

tmp.cor.melt.df <- filter(RFD.cor.melt.df, clone %in% c("410",mt.ids))
tmp.spear.df <- filter(RFD.spear.df, clone %in% c("410",mt.ids))

max.part <- 0.7                
supp.fig.partition.scatter.plt <- tmp.cor.melt.df %>%
  mutate(clone = factor(clone, levels = names(clone.labs))) %>%
  ggplot(aes(x = RFD, y = Partition)) +
  geom_hex(bins=120) + ggtitle(mk) +
  scale_fill_gradientn(colours = (brewer.pal(n=9,name="Blues")[2:8])) +
  geom_vline(xintercept = 0, colour = "grey70", size = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", size = 0.5) +
  scale_y_continuous(breaks=c(-0.6, 0, 0.6), 
                     limits = c(-max.part,max.part)) +
  scale_x_continuous(breaks=c(-0.6, 0, 0.6), 
                     limits =  c(-max.part,max.part)) +
  geom_smooth(se=F,method="lm",size=0.3,colour="red") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  geom_text(data = tmp.spear.df, 
            aes(x = -0.4, y = 0.5, 
                label = tex.label),
            parse = TRUE, size = 5) +
  facet_grid(clone~timepoint, labeller = labeller(.rows = clone.labs,
                                                  .cols = tp.labs)) +
  theme(
    panel.spacing = unit(0.3, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 7 * .pt),
    axis.text = element_text(size = 6 * .pt),
    legend.text = element_text(size = 6 * .pt),
    strip.text.x = element_text(size = 5 * .pt),
    strip.text.y = element_text(size = 5 * .pt),
    strip.background = element_rect(colour="black", fill="#e3e3e3"),
    aspect.ratio = 0.8)

ggsave(filename = file.path(plots.dir, paste0("FigS4I_SCAR_Inputs_scatter.pdf")),
       plot =   supp.fig.partition.scatter.plt,
       width = 16, height = 12, device = cairo_pdf)


load(file.path(scar.data.dir,"SCAR_H3K27ac_bins_stranded_w250.RData"))

all.csaw.bins.df$orient_IZ <- "Upstream"
all.csaw.bins.df$orient_IZ[all.csaw.bins.df$dist.to.IZ >= 0] <- "Downstream"

all.csaw.bins.df$strand_IZ <- "Lagging"

all.csaw.bins.df$strand_IZ[(all.csaw.bins.df$strand == "+") & 
                             (all.csaw.bins.df$orient_IZ == "Downstream")] <- "Leading"

all.csaw.bins.df$strand_IZ[(all.csaw.bins.df$strand == "-") & 
                             (all.csaw.bins.df$orient_IZ == "Upstream")] <- "Leading"


all.csaw.bins.df$strand_IZ[(all.csaw.bins.df$strand == "-") & 
                             (all.csaw.bins.df$orient_IZ == "Downstream")] <- "Lagging"

all.csaw.bins.df$strand_IZ[(all.csaw.bins.df$strand == "-") & 
                             (all.csaw.bins.df$orient_IZ == "Upstream")] <- "Lagging"


all.csaw.bins.df$peak.group <- "None"
all.csaw.bins.df$peak.group[all.csaw.bins.df$H3K27ac] <- "H3K27ac"
all.csaw.bins.df$peak.group[all.csaw.bins.df$H3K27me3] <- "H3K27me3"

all.csaw.bins.df$peak.group <- factor(all.csaw.bins.df$peak.group, levels = c("H3K27ac","H3K27me3","None"))
all.csaw.bins.df$condition <- clone.conds[as.character(all.csaw.bins.df$clone)]



all.csaw.wide.df <- c()
for (cid in c("410","439","442")) {
  for (tid in c("T0","T1","T3","T8")) {
    tmp.bins.df <- filter(all.csaw.bins.df, clone == cid, timepoint == tid, is.good)
    tmp.annot.df <- unique(tmp.bins.df[,c("bins","peak.group", "orient_IZ",
                                          "clone","is.good",
                                          "timepoint","condition")])
    
    tmp.bins.mat <- reshape2::dcast(data = tmp.bins.df, bins ~ strand, fill = 0, value.var = "CPM")
    colnames(tmp.bins.mat) <- c("bins","F.cpm", "R.cpm")
    
    tmp.annot.df <- full_join(tmp.annot.df, tmp.bins.mat)
    
    tmp.bins.mat <- reshape2::dcast(data = tmp.bins.df, bins ~ strand, fill = 0, value.var = "log2FC")
    colnames(tmp.bins.mat) <- c("bins","F.log2FC", "R.log2FC")
    
    tmp.annot.df <- full_join(tmp.annot.df, tmp.bins.mat)
    
    
    tmp.annot.df$Partition <- (tmp.annot.df$F.cpm - tmp.annot.df$R.cpm) /
      (tmp.annot.df$F.cpm + tmp.annot.df$R.cpm)
    
    tmp.annot.df$Partition.FC <- log2((tmp.annot.df$F.cpm + 0.01) / (tmp.annot.df$R.cpm + 0.01)) 
    
    
    tmp.annot.df$orient_IZ <- factor(tmp.annot.df$orient_IZ, levels = c("Upstream","Downstream"))
    all.csaw.wide.df <- rbind(all.csaw.wide.df, tmp.annot.df)
  }
}

CPM.cutoff <- 0.3
all.csaw.wide.df$Partition.adj <- ifelse(all.csaw.wide.df$orient_IZ == "Upstream",
                                         - all.csaw.wide.df$Partition,
                                         all.csaw.wide.df$Partition)

lb.bins <- filter(all.csaw.wide.df,(F.cpm + R.cpm) > CPM.cutoff) %>% 
  group_by(clone, timepoint, peak.group,orient_IZ,is.good) %>%
  tally() %>% as.data.frame()

all.csaw.wide.df$condition <- factor(all.csaw.wide.df$condition, 
                                     levels = c("WT","MCM2-2A"))


clone.lab.colors <- c("WT" = "#0C46A0FF",
                  "MCM2-2A #1" = "#FFA626FF",
                  "MCM2-2A #2" = "darkorange")

all.csaw.wide.df$clone.labs <-
  factor(clone.labs[as.character(all.csaw.wide.df$clone)],
         levels = c("WT","MCM2-2A #1","MCM2-2A #2"))




scar.fig.S4F.plt <- filter(all.csaw.wide.df, (F.cpm + R.cpm) > CPM.cutoff,
                      is.good) %>%
  ggplot(aes(x = peak.group, y = Partition.adj, 
             fill = clone.labs)) + 
  geom_boxplot(outlier.shape = NA,
               notch = 1, lwd = 0.8) + 
  
  scale_fill_manual(values=clone.lab.colors) +
  scale_y_continuous(breaks = c(-0.6, 0, 0.6), limits = c(-1.1,1.1)) +
  ylab("H3K27ac Partition") + xlab("Genomic Regions") +
  theme_bw(base_size = 20, base_family = "Helvetica") + 
  facet_wrap(~ timepoint,nrow = 1,labeller = labeller(.cols = tp.labs)) +
  labs(fill = "", colour = "Strand", alpha = "Strand") +
  guides(alpha = guide_legend(title = "Relative to IZ",
                              override.aes = list(fill = c("Upstream"=part.cond.cols["MCM2-2A"],
                                                           "Downstream"=part.cond.cols["MCM2-2A"])))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 7 * .pt),
        axis.text = element_text(size = 5 * .pt),
        legend.text = element_text(size = 6 * .pt),
        panel.spacing = unit(0.35, "lines"),
        strip.text.x = element_text(size = 5 * .pt),
        strip.text.y = element_text(size = 5 * .pt),
        strip.background = element_rect(colour="black", fill="#e3e3e3"))

ggsave(filename = file.path(plots.dir, paste0("FigS4F_SCAR_H3K27ac_partition_peaks.pdf")),
       plot = scar.fig.S4F.plt,
       width = 16, height = 14, device = cairo_pdf)


tmp1.df <- filter(all.csaw.wide.df, (F.cpm + R.cpm) > 0.3,
       is.good) %>% mutate(grouping.time =
                             factor(paste0(peak.group,"_",timepoint))) %>%
  mutate(clone = factor(clone, levels = c("410","439","442"))) %>%
  as.data.frame()


tmp2.df <- filter(all.csaw.wide.df, (F.cpm + R.cpm) > 0.3,
                  is.good) %>% mutate(grouping.clone =
                                        factor(paste0(peak.group,"_",clone))) %>%
  mutate(timepoint = factor(timepoint, levels = c("T0","T1","T3","T8"))) %>%
  as.data.frame()


tmp1.test.df <- full_join(compare_means(formula = Partition.adj ~ clone,
              data = tmp1.df, group.by = "grouping.time",
              method = "wilcox.test") %>%
  mutate(peak.group = strsplit2(as.character(grouping.time),"_")[,1],
         timepoint = strsplit2(as.character(grouping.time),"_")[,2]) %>%
    as.data.frame(), tmp1.df %>% group_by(grouping.time) %>%
  wilcox_effsize(Partition.adj ~ clone) %>%
  as.data.frame())

tmp1.test.df <- tmp1.test.df[,c("group1","group2",
                                "n1","n2","peak.group","timepoint",
                                "p","p.adj","p.format",
                                "p.signif","method",
                                "effsize","magnitude")]

  
tmp2.test.df <- full_join(compare_means(formula = Partition.adj ~ timepoint,
                              data = tmp2.df, group.by = "grouping.clone",
                              method = "wilcox.test") %>%
  mutate(peak.group = strsplit2(as.character(grouping.clone),"_")[,1],
         clone = strsplit2(as.character(grouping.clone),"_")[,2]) %>%
  as.data.frame(),tmp2.df %>% group_by(grouping.clone) %>%
    wilcox_effsize(Partition.adj ~ timepoint) %>%
    as.data.frame())


tmp2.test.df <- tmp2.test.df[,c("group1","group2",
                                "n1","n2","peak.group","clone",
                                "p","p.adj","p.format",
                                "p.signif","method",
                                "effsize","magnitude")]

write.xlsx(list("Clone" = tmp1.test.df,
                "Timepoint" = tmp2.test.df),
           file = file.path(scar.data.dir, paste0("FigS4F_SCAR_H3K27ac_partition_peaks.xlsx")))




