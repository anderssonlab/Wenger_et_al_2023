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
scar.data.dir <- file.path(data.dir,"SCAR_seq")

external.data.dir <- "data/external"

## Load process external data and functions
load(file.path(external.data.dir, "Processed_external.RData"))
source(file.path(wd.dir,"utils/utils.R"))
num.cores <- detectCores() / 2

chroms <- paste("chr",1:19,sep="") 

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

load(file.path(chip.data.dir,"ChIP_processed_peaks.RData"))

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
scar.peaks.list[["ChIP_H2AK119ub1_strict"]] <- chip.external.peaks.list$H2AK119ub1[chip.external.peaks.list$H2AK119ub1$pval >= -log10(0.05)]



Okazaki.df <- as.data.frame(OK.gr)
Okazaki.df$ID <- "OK-seq"
Okazaki.df$sample <- "OK-seq"
Okazaki.df$mark <- NA
Okazaki.df$condition <- NA
Okazaki.df$clone <- "263"
Okazaki.df$timepoint <- NA
Okazaki.df$replicate <- NA
Okazaki.df$project <- NA

rfd.dir.list <- list("MCM2" = "data/SCAR_seq/MCM2/rfd_files",
                     "Crosslinked" = "data/SCAR_seq/Crosslinked/rfd_files",
                     "Rescue" = "data/SCAR_seq/Rescue/rfd_files",
                     "Pole4" = "data/SCAR_seq/Pole4/rfd_files",
                     "Input" = "data/SCAR_seq/Input/rfd_files",
                     "InputPooled" = "data/SCAR_seq/InputPooled/rfd_files")

clone.condition <- c("410" = "WT", "439" = "MCM2-2A","442" = "MCM2-2A", 
                     "438" = "MCM2-2A",
                     "423" = "MCM2-Y90A",
                     "588" = "MCM2-2A-R",
                     "551" = "Pole4-2A", "552" = "Pole4-2A")

full.design.df <- read.csv(file = file.path(scar.data.dir,"RFD_files_info.tsv"),
                           header = TRUE, sep = "\t")



full.design.df$ID <- full.design.df$sample
full.design.df$clone <- factor(full.design.df$clone, levels = names(clone.conds))
full.design.df$mark <- factor(full.design.df$mark)
full.design.df$replicate <- factor(full.design.df$replicate)
full.design.df$timepoint <- factor(full.design.df$timepoint)
full.design.df$project <- factor(full.design.df$project)



for (prj in setdiff(unique(full.design.df$project), c("MCM2","Crosslinked"))) {
  tmp.design.df <- filter(full.design.df, project == prj)
  for (mk in unique(gsub("X","",tmp.design.df$mark))) {
    
    run.id <- paste0(prj, "_",mk)
    print(run.id)
    
    scar.design.df <- filter(tmp.design.df, mark %in% c(mk,"inputX"))
    rownames(scar.design.df) <- scar.design.df$sample
    
    
    cls <- c("seqnames","start","end","F","R","F.cpm","R.cpm","RFD.raw","RFD",
             "RFD.deriv","score","zero.deriv")
    cls.select <- c(1:5,13,14,8,15,16,17,18)
    if (prj %in% c("Input")) {
      cls.select <- 1:12
    }
    
    print("Loading partition files...")
    cl <- makeCluster(min(nrow(scar.design.df),30), type = "FORK")
    registerDoParallel(cl)
    
    SCAR.full.df <- foreach (i = 1:nrow(scar.design.df), .combine = "rbind",
                             .packages = c("GenomicRanges", "reshape2",
                                           "dplyr","stringr", "stringi")) %dopar% {
                                             
                                             SCAR.df <- read.table(scar.design.df[i,]$file,
                                                                   sep="\t",header=FALSE,
                                                                   as.is=TRUE, fill = TRUE)[,cls.select]
                                             
                                             colnames(SCAR.df) <- cls
                                             SCAR.df <- dplyr::filter(SCAR.df, seqnames %in% chroms)
                                             SCAR.df$names <- paste(SCAR.df$seqnames,SCAR.df$start,sep=":") %>% paste(.,SCAR.df$end,sep="-") # genomic identifier
                                             SCAR.df$exprs <- SCAR.df[,"F"] + SCAR.df[,"R"] # total raw counts
                                             SCAR.df$CPM <- SCAR.df$F.cpm + SCAR.df$R.cpm # total counts pr. million
                                             
                                             SCAR.gr <- makeGRangesFromDataFrame(SCAR.df, keep.extra.columns = TRUE)
                                             bl.ovl <- findOverlaps(SCAR.gr, mm10.blacklist.noovrl.gr, minoverlap = 1,
                                                                    select = "all", type = "any")
                                             SCAR.gr <- SCAR.gr[-unique(queryHits(bl.ovl))]
                                             SCAR.gr$IZ <- NA
                                             
                                             tmp.df <- as.data.frame(SCAR.gr)
                                             for (fts in setdiff(colnames(scar.design.df),"file")) {
                                               tmp.df[,fts] <- scar.design.df[i,fts]
                                             }
                                             return(tmp.df)
                                           }
    stopCluster(cl)
    
    
    print("Annotating regions...")
    RFD.gr <- makeGRangesFromDataFrame(rbind(Okazaki.df, SCAR.full.df),
                                       keep.extra.columns = TRUE)
    
    dist.RT <- as.data.frame(distanceToNearest(RFD.gr, RT.gr))
    RFD.gr$RT <- RT.gr[dist.RT$subjectHits]$RT.class
    if (FALSE) {
      for (pkn in names(scar.peaks.list)) {#names(scar.peaks.list)) {
        print(pkn)
        tmp.gr <- scar.peaks.list[[pkn]]
        tmp.gr <- tmp.gr[seqnames(tmp.gr) %in% chroms]
        seqlevels(tmp.gr) <- chroms
        tmp.gr <- keepSeqlevels(tmp.gr, chroms)
        rfd.dist.chip <- as.data.frame(distanceToNearest(resize(RFD.gr, width = 2, fix = "center"), 
                                                         tmp.gr))
        
        elementMetadata(RFD.gr)[[pkn]] <- FALSE
        elementMetadata(RFD.gr)[[paste0(pkn,"_ovrl")]] <- FALSE
        elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance <= 100000]])[[pkn]] <- TRUE
        elementMetadata(RFD.gr[rfd.dist.chip$queryHits[rfd.dist.chip$distance == 0]])[[paste0(pkn,"_ovrl")]] <- TRUE
      }
    }
    
    
    rm(list = c("SCAR.full.df", "rfd.dist.chip"))
    gc()
    
    CPM.cutoff <- 0.3
    print("Computing correlations with OK-seq...")
    RFD.cor.mat <- as.data.frame(RFD.gr) %>% 
      filter(mark == mk | sample == "OK-seq", (F.cpm + R.cpm) >= CPM.cutoff) %>%
      group_by(names, clone, mark, timepoint, condition) %>% 
      summarise(RFD = mean(RFD, na.rm = TRUE)) %>% as.data.frame()
    
    RFD.cor.mat$ID <- paste(RFD.cor.mat$clone, RFD.cor.mat$mark, RFD.cor.mat$timepoint, sep = "_")
    RFD.cor.mat[RFD.cor.mat$clone == "263",]$ID <- "OK-seq"
    
    RFD.cor.df  <- reshape2::dcast(RFD.cor.mat, names ~ ID, value.var = "RFD", fun.aggregate = mean)
    
    cor.comp <- colnames(RFD.cor.df)[-1]
    cor.comp <- cor.comp[order(cor.comp)]
    #cor.comp <- c("OK_WT",cor.comp)
    cor.df <- c()
    for (i in 1:(length(cor.comp)-1)) {
      for (j in (i+1):length(cor.comp)) {
        comp.x <- cor.comp[i]
        comp.y <- cor.comp[j]
        cor.obj <- cor.test(RFD.cor.df[,comp.x], RFD.cor.df[,comp.y], method = "spearman",
                            use = "pairwise.complete.obs")
        tmp.df <- data.frame(cor.x = comp.x, cor.y = comp.y, 
                             spearman.cor = cor.obj$estimate,
                             cor.pvalue = cor.obj$p.value,
                             stringsAsFactors = FALSE)
        cor.df <- rbind(cor.df, tmp.df)
      }
    }
    
    RFD.spear.df <- filter(cor.df, cor.y == "OK-seq")
    RFD.spear.df <- RFD.spear.df[,-2]
    
    colnames(RFD.spear.df) <- c("SCAR","spearman.cor","cor.pvalue")
    
    RFD.spear.df$timepoint <- factor(strsplit2(RFD.spear.df$SCAR, "_")[,3])
    
    RFD.spear.df$tex.label <- paste0("$\\rho = ", sprintf("%0.4f",RFD.spear.df$spearman.cor), "$")
    RFD.spear.df$tex.label <- paste0("rho == ", sprintf("%0.4f",RFD.spear.df$spearman.cor))
    
    RFD.cor.melt.df <- melt(RFD.cor.df, id.vars = c("names","OK-seq"),
                            variable.name = "SCAR", value.name = "RFD") %>%
      filter(!is.nan(RFD) & !is.nan('OK-seq'))
    
    
    colnames(RFD.cor.melt.df) <- c("names", "RFD", "SCAR", "Partition")
    
    RFD.cor.melt.df$timepoint <-  factor(strsplit2(RFD.cor.melt.df$SCAR, "_")[,3])
    
    fout <- paste0("SCAR_rfd_" ,run.id, ".RData")
    save(file = file.path(scar.data.dir,fout),
         list = c("RFD.gr","RFD.cor.melt.df","RFD.spear.df"))
  }
}

#### For H3K27ac binned analyses

### Analysis of K27ac signal over peaks

outfile <- paste0("SCAR_rfd_MCM2", "_","H3K27ac",".RData")
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


overlap.pairs <- findOverlaps(OK.ext.gr, RFD.gr)
RFD.break.all.gr <- RFD.gr[subjectHits(overlap.pairs)]
RFD.break.all.gr$break_ID <- OK.ext.gr$names[queryHits(overlap.pairs)]
RFD.break.all.gr$RT <- OK.ext.gr$RT[queryHits(overlap.pairs)]
RFD.break.all.gr$dist <- start(RFD.break.all.gr) - OK.ext.gr$break_start[queryHits(overlap.pairs)]

RFD.break.ok.df <- RFD.break.all.gr %>% as.data.frame() %>% 
  filter(sample == "OK-seq")


RFD.break.all.gr$H3K27ac <- FALSE
RFD.break.all.gr$H3K27me3 <- FALSE

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

OK.ext.gr <- sortSeqlevels(OK.ext.gr)
OK.ext.gr <- sort(OK.ext.gr)

names(RFD_max.gr) <- RFD_max.gr$break_ID
names(RFD_min.gr) <- RFD_min.gr$break_ID

RFD_max.gr <- RFD_max.gr[OK.ext.gr$names]
RFD_min.gr <- RFD_min.gr[OK.ext.gr$names]


scar.bam.dir <- file.path(scar.data.dir, "MCM2/bam_files")
scar.stats.df <- read.csv(file.path(scar.data.dir, "SCAR_stats.tsv"),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(scar.stats.df) <- scar.stats.df$ID
scar.stats.df$total.endo <- scar.stats.df$endo.plus.reads + scar.stats.df$endo.minus.reads
scar.stats.df$total.exo <- scar.stats.df$exo.plus.reads + scar.stats.df$exo.minus.reads


tmp.bams <- list.files(scar.bam.dir, pattern = "nodup.mm10.bam$")

scar.bams <- tmp.bams[grepl("H3K27ac_",tmp.bams)]
names(scar.bams) <- strsplit2(scar.bams,"\\.")[,1]

input.bams <- tmp.bams[grepl("H3K27acinput_",tmp.bams)]
names(input.bams) <- strsplit2(input.bams,"\\.")[,1]

all.csaw.bins.df <- c()
for(cid in c("410","439","442")) {
  k27ac.gr <- chip.peaks.list[[paste0("H3K27ac_",cid)]]
  k27me3.gr <- chip.peaks.list[[paste0("H3K27me3_",cid)]]
  
  k27ac.excl.gr <- k27ac.gr[!overlapsAny(k27ac.gr,k27me3.gr, minoverlap = 1)]
  k27me3.excl.gr <- k27me3.gr[!overlapsAny(k27me3.gr,k27ac.gr, minoverlap = 1)]
  
  for(tp in c("T0","T1","T3","T8")) {
    print(paste(cid, tp))
    csaw.bams <- c(scar.bams[grepl(cid,scar.bams) & grepl(tp,scar.bams)],
                   input.bams[grepl(cid,input.bams) & grepl(tp,input.bams)])
    
    lib.sizes <- scar.stats.df[names(csaw.bams),]$total.endo
    
    
    read.param <- readParam(discard = mm10.blacklist.noovrl.gr,
                            restrict = setdiff(chroms,c("chrX","chrY")), 
                            pe = "first", forward = NULL)
    
    
    csaw.chip.bkg <- strandedCounts(file.path(scar.bam.dir,csaw.bams), 
                                    param = read.param, 
                                    BPPARAM = MulticoreParam(workers = length(csaw.bams)),
                                    width = 10000, 
                                    bin = TRUE, ext = 150, shift = 75)
    
    
    csaw.chip.bin <- strandedCounts(file.path(scar.bam.dir,csaw.bams), 
                                    param = read.param, 
                                    BPPARAM = MulticoreParam(workers = 32),
                                    width = 250, 
                                    bin = TRUE, ext = 150, shift = 75)
    
    colnames(csaw.chip.bin) <- names(csaw.bams)
    
    scale.info <- scaleControlFilter(csaw.chip.bkg[,1:2], 
                                     csaw.chip.bkg[,3:4], 
                                     assay.data = 1,assay.back = 1)
    
    filtered.unspiked <- filterWindowsControl(csaw.chip.bin[,1:2], 
                                              csaw.chip.bin[,3:4], 
                                              assay.data = 1,assay.back = 1,
                                              scale.info = scale.info) 
    
    rownames(csaw.chip.bin) <- paste0(seqnames(rowRanges(csaw.chip.bin)),":",
                                      start(rowRanges(csaw.chip.bin)), "-",
                                      end(rowRanges(csaw.chip.bin)), ",",
                                      strand(rowRanges(csaw.chip.bin)))
    
    
    keep.unspiked <- filtered.unspiked$filter > log2(1.5)
    
    csaw.chip.dge <- asDGEList(csaw.chip.bin, assay.id = 1, lib.size = lib.sizes)
    
    csaw.cpms <- cpm(csaw.chip.dge, log = FALSE, lib.size = lib.sizes)
    colnames(csaw.cpms) <- colnames(csaw.chip.bin)
    
    csaw.chip.gr <- rowRanges(csaw.chip.bin)
    
    dist.IZ.df <- as.data.frame(distanceToNearest(resize(csaw.chip.gr, width = 1, fix = "center"),
                                                  resize(OK.ext.gr, width = 1, fix = "center")))
    
    dist.max.df <- as.data.frame(distanceToNearest(resize(csaw.chip.gr, width = 1, fix = "center"),
                                                   resize(RFD_max.gr, width = 1, fix = "center")))
    
    dist.min.df <- as.data.frame(distanceToNearest(resize(csaw.chip.gr, width = 1, fix = "center"),
                                                   resize(RFD_min.gr, width = 1, fix = "center")))
    
    dist.IZ.max <- start(resize(RFD_max.gr, fix = "center", width = 1)) -
      start(resize(OK.ext.gr, fix = "center", width = 1))
    
    dist.IZ.min <- start(resize(OK.ext.gr, fix = "center", width = 1)) -
      start(resize(RFD_min.gr, fix = "center", width = 1))
    
    dist.to.IZ <- start(resize(csaw.chip.gr,fix = "center", width = 1)[dist.IZ.df$queryHits]) -
      start(resize(OK.ext.gr,fix = "center", width = 1)[dist.IZ.df$subjectHits])
    
    dist.to.MAX <- start(resize(csaw.chip.gr,fix = "center", width = 1)[dist.max.df$queryHits]) -
      start(resize(RFD_max.gr,fix = "center", width = 1)[dist.max.df$subjectHits])
    
    dist.to.MIN <- start(resize(csaw.chip.gr,fix = "center", width = 1)[dist.min.df$queryHits]) -
      start(resize(RFD_min.gr,fix = "center", width = 1)[dist.min.df$subjectHits])
    
    
    good.zone <- (dist.IZ.df$distance >= 10000) & ((dist.max.df$distance <= dist.IZ.max[dist.max.df$subjectHits]) |
                                                     (dist.min.df$distance <= dist.IZ.min[dist.min.df$subjectHits]))
    
    
    csaw.bins.df <- data.frame(bins = paste0(seqnames(csaw.chip.gr),":",
                                             start(csaw.chip.gr),"-",
                                             end(csaw.chip.gr)),
                               strand = strand(csaw.chip.gr),
                               log2FC = filtered.unspiked$filter,
                               CPM.raw = rowMeans(csaw.cpms[,1:2]),
                               CPM = rowMeans(pmax((csaw.cpms[,1:2] - csaw.cpms[,3:4]), 0)),
                               H3K27ac = overlapsAny(csaw.chip.gr,k27ac.excl.gr, 
                                                     minoverlap = 125),
                               H3K27me3 = overlapsAny(csaw.chip.gr,k27me3.excl.gr, 
                                                      minoverlap = 125),
                               dist.to.IZ = dist.to.IZ,
                               dist.to.MAX = dist.to.MAX,
                               dist.to.MIN = dist.to.MIN,
                               is.good = good.zone,
                               clone = cid,
                               timepoint = tp,
                               stringsAsFactors = FALSE)
    
    csaw.bins.df <- csaw.bins.df[csaw.bins.df$log2FC > 0,]
    
    all.csaw.bins.df <- rbind(all.csaw.bins.df, csaw.bins.df)
    
  }
}

save(file = file.path(scar.data.dir,"SCAR_H3K27ac_bins_stranded_w250.RData"),
     list = c("all.csaw.bins.df","scar.stats.df"))






