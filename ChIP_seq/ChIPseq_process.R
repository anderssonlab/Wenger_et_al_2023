library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)
library(stringr)
library(ChIPseeker)
library(stringi)
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
library(edgeR)

## ======================= Load annotations ====================================
wd.dir <- "Wenger_et_al_2022"
setwd(wd.dir)

external.data.dir <- "data"

## Load process external data and functions
load(file.path(wd.dir,external.data.dir,"external", "Processed_external.RData"))
source(file.path("utils","utils.R"))

num.cores <- detectCores() / 2

register(MulticoreParam(workers = num.cores))

data.dir <- "data"
external.data.dir <- file.path(data.dir,"external")
annot.version <- "gencode.vM23"
assembly <- "mm10"

chrom.sizes.file <- file.path(data.dir, paste0(assembly, ".chrom.sizes"))
chrom.sizes.df <- read.csv2(file = chrom.sizes.file, header = FALSE, col.names = c("chrom", "length"), sep = "\t")

annot.data.gtf <- file.path(external.data.dir, paste0(annot.version,".annotation.gtf.gz"))
annot.sqlite.file <- gsub(".gtf.gz",".sqlite", annot.data.gtf)

if (!file.exists(annot.sqlite.file)) {
  gencode.txdb <- makeTxDbFromGFF(file = annot.data.gtf, format="gtf", chrominfo=chrom.sizes.df, 
                                  dataSource=paste("arke:", annot.data.gtf, sep = ""),
                                  organism="Mus musculus")
  saveDb(gencode.txdb, file = annot.sqlite.file)
} else {
  gencode.txdb <- loadDb(annot.sqlite.file)
}


## ===================  Load external data =====================================
chroms <- c(paste("chr",1:19,sep=""), "chrX", "chrY")


### Load file tsv specifying path to bam files
chip.external.data.dir <- "data/ChIP_seq"
chip.stats.df <- read.csv(file.path(chip.external.data.dir,"chip_stats.tsv"),
                          sep = "\t", header = TRUE)

rownames(chip.stats.df) <- chip.stats.df$ID

chip.external.stats.df <- read.csv(file.path(chip.external.data.dir,"chip_external_stats.tsv"),
                                   sep = "\t", header = TRUE)

rownames(chip.external.stats.df) <- chip.external.stats.df$ID

chip.rescue.stats.df <- read.csv(file.path(chip.external.data.dir,"chip_rescue_stats.tsv"),
                                   sep = "\t", header = TRUE)

rownames(chip.rescue.stats.df) <- chip.rescue.stats.df$ID


#### Read in raw peaks
chip.data.dir <- file.path(data.dir, "ChIP_seq")

cls <- c("seqnames", "start", "end", "name", "score", "strand",
         "log2FC", "pval", "qval", "summit")


chip.broad.dir <- file.path(chip.data.dir,"MAPPING_mm10/peak_files/broad")
chip.broad.files <- list.files(chip.broad.dir, pattern = "broadPeak.gz")

names(chip.broad.files) <- strsplit2(chip.broad.files, "\\.")[,1]

chip.broad.raw.peaks.list <- mclapply(names(chip.broad.files), 
                                mc.cores = min(num.cores,length(chip.broad.files)),
                                FUN = function(pkf) {
  pk.df <- read.csv(file.path(chip.broad.dir,chip.broad.files[pkf]), 
                    sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE)
  colnames(pk.df) <- cls[1:9]
  pk.df <- pk.df[pk.df$seqnames %in% chroms,]
  pk.gr <- makeGRangesFromDataFrame(pk.df, keep.extra.columns = TRUE)
  pk.gr <- pk.gr[!overlapsAny(pk.gr, mm10.blacklist.noovrl.gr)]
  return(pk.gr)
})

names(chip.broad.raw.peaks.list) <- names(chip.broad.files)

chip.narrow.dir <- file.path(chip.data.dir,"MAPPING_mm10/peak_files/narrow")
chip.narrow.files <- list.files(chip.narrow.dir, pattern = "narrowPeak.gz")
names(chip.narrow.files) <- strsplit2(chip.narrow.files, "\\.")[,1]


chip.narrow.raw.peaks.list <- mclapply(names(chip.narrow.files), 
                                      mc.cores = min(num.cores,length(chip.narrow.files)),
                                      FUN = function(pkf) {
                                        pk.df <- read.csv(file.path(chip.narrow.dir,chip.narrow.files[pkf]), 
                                                          sep = "\t", header = FALSE,
                                                          stringsAsFactors = FALSE)
                                        colnames(pk.df) <- cls
                                        pk.df <- pk.df[pk.df$seqnames %in% chroms,]
                                        pk.gr <- makeGRangesFromDataFrame(pk.df, keep.extra.columns = TRUE)
                                        pk.gr <- pk.gr[!overlapsAny(pk.gr, mm10.blacklist.noovrl.gr)]
                                        return(pk.gr)
                                      })


names(chip.narrow.raw.peaks.list) <- names(chip.narrow.files)
chip.raw.peaks.list <- c(chip.broad.raw.peaks.list, chip.narrow.raw.peaks.list)

rm(list = c("chip.broad.raw.peaks.list","chip.narrow.raw.peaks.list"))


### Create consensus peaks from biological replicates

chip.peaks.list <- list()
#chip.reprod2.peaks.list <- list()
for (mk in unique(chip.stats.df$mark)) {
  pk.mk.list <- chip.raw.peaks.list[grepl(mk,names(chip.raw.peaks.list))]
  for (cid in unique(chip.stats.df$clone)) { #unique(chip.stats.df$clone)) {
    rid <- paste0(mk, "_",cid)
    print(rid)
    if (sum(grepl(cid,names(pk.mk.list))) == 0) {
      next
    }
    pk.clone.list <- pk.mk.list[grepl(cid,names(pk.mk.list))]
    pk.pr.gr <- pk.clone.list[grepl("_pr",names(pk.clone.list))][[1]]
    pk.rep.list <- pk.clone.list[!grepl("_pr",names(pk.clone.list))]
    tmp.gr <- GRanges()
    for(i in 1:length(pk.rep.list)) {
      tmp.gr <- c(tmp.gr, pk.rep.list[[i]])  
    }
    
    pk.pooled.gr <- GenomicRanges::reduce(tmp.gr)
    pr.mat <- matrix(FALSE, nrow = length(pk.pr.gr), ncol = length(pk.rep.list))
    pooled.mat <- matrix(FALSE, nrow = length(pk.pooled.gr), ncol = length(pk.rep.list))
    for (j in 1:length(pk.rep.list)) {
      pr.mat[,j] <- overlapsAny(pk.pr.gr, pk.rep.list[[j]], minoverlap = 1)
      pooled.mat[,j] <- overlapsAny(pk.pooled.gr, pk.rep.list[[j]], minoverlap = 1)
    }
    
    if (mk == "SUZ12") {
      chip.peaks.list[[rid]] <-  annotateRegions(pk.pr.gr[rowSums(pr.mat) >= 2], gencode.txdb, 
                                                 tss.region = c(-3000,3000),
                                                 Enh.gr, SEnh.gr, cpg.islands.gr,
                                                 all_CpG_sites, cpg.meth.gr)
      
    } else {
      chip.peaks.list[[rid]] <- annotateRegions(pk.pooled.gr[rowSums(pooled.mat) >= 2], gencode.txdb, 
                                                        tss.region = c(-3000,3000),
                                                        Enh.gr, SEnh.gr, cpg.islands.gr,
                                                        all_CpG_sites, cpg.meth.gr)
    }
  }
}


chip.rescue.peaks.dir <- "data/ChIP_seq/MAPPING_mm10/peak_files/broad/rescue"

chip.rescue.broad.files <- list.files(chip.rescue.peaks.dir, pattern = "broadPeak.gz")

names(chip.rescue.broad.files) <- strsplit2(chip.rescue.broad.files, "\\.")[,1]

chip.raw.rescue.peaks.list <- mclapply(names(chip.rescue.broad.files), 
         mc.cores = min(num.cores,length(chip.rescue.broad.files)),
         FUN = function(pkf) {

  pk.df <- read.csv(file.path(chip.rescue.peaks.dir,chip.rescue.broad.files[pkf]), sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE)
  colnames(pk.df) <- cls[1:9]
  pk.df <- pk.df[pk.df$seqnames %in% chroms,]
  pk.gr <- makeGRangesFromDataFrame(pk.df, keep.extra.columns = TRUE)
  pk.gr <- pk.gr[!overlapsAny(pk.gr, mm10.blacklist.noovrl.gr)]
  return(pk.gr)
})

names(chip.raw.rescue.peaks.list) <- names(chip.rescue.broad.files)

chip.rescue.peaks.list <- list()
for (mk in unique(chip.rescue.stats.df$mark)) {
  pk.mk.list <- chip.raw.rescue.peaks.list[grepl(mk,names(chip.raw.rescue.peaks.list))]
  for (cid in unique(chip.rescue.stats.df$clone)) { #unique(chip.stats.df$clone)) {
    rid <- paste0(mk, "_",cid)
    print(rid)
    if (sum(grepl(cid,names(pk.mk.list))) == 0) {
      next
    }
    pk.clone.list <- pk.mk.list[grepl(cid,names(pk.mk.list))]
    pk.pr.gr <- pk.clone.list[grepl("_pr",names(pk.clone.list))][[1]]
    pk.rep.list <- pk.clone.list[!grepl("_pr",names(pk.clone.list))]
    tmp.gr <- GRanges()
    for(i in 1:length(pk.rep.list)) {
      tmp.gr <- c(tmp.gr, pk.rep.list[[i]])  
    }
    
    pk.pooled.gr <- GenomicRanges::reduce(tmp.gr)
    pr.mat <- matrix(FALSE, nrow = length(pk.pr.gr), ncol = length(pk.rep.list))
    pooled.mat <- matrix(FALSE, nrow = length(pk.pooled.gr), ncol = length(pk.rep.list))
    for (j in 1:length(pk.rep.list)) {
      pr.mat[,j] <- overlapsAny(pk.pr.gr, pk.rep.list[[j]], minoverlap = 1)
      pooled.mat[,j] <- overlapsAny(pk.pooled.gr, pk.rep.list[[j]], minoverlap = 1)
    }
    
  
    chip.rescue.peaks.list[[rid]] <- annotateRegions(pk.pooled.gr[rowSums(pooled.mat) >= 2], gencode.txdb, 
                                                      tss.region = c(-3000,3000),
                                                      Enh.gr, SEnh.gr, cpg.islands.gr,
                                                      all_CpG_sites, cpg.meth.gr)
  }
}


### Load external chip files


chip.external.peak.dir <- "data/ChIP_seq/external/peak_files"
chip.ext.files <- c("H2AK119ub1" = "ChIP_Blackledge_H2AK119ub1_pr_R1.peaks.broadPeak.gz",
                    "H3K36me3" = "ChIP_EncodePennstate_H3K36me3_pr.peaks.broadPeak.gz",
                    "H3K27me1" = "ChIP_Healy_H3K27me1_r1.peaks.broadPeak.gz",
                    "H3K27me2"= "ChIP_Healy_H3K27me2_r1.peaks.broadPeak.gz",
                    "H3K9me2" = "ChIP_Jiang_H3K9me2_J_r1.peaks.broadPeak.gz",
                    "H3K36me2" = "ChIP_Tuberfield_H3K36me2_pr_R1.peaks.broadPeak.gz",
                    "H3K4me1" = "ChIP_Yan_H3K4me1_pr.peaks.broadPeak.gz",
                    "SUV39H1" = "ChIP_Bulut_SUV39H1_pr_R1.peaks.broadPeak.gz",
                    "SUV39H2" = "ChIP_Bulut_SUV39H2_pr_R1.peaks.broadPeak.gz",
                    "SETDB1" = "ChIP_Bilodeu_SETDB1_pr.peaks.broadPeak.gz")


chip.external.peaks.list <- list()
cls2 <- cls[1:9] 
for (mk in names(chip.ext.files)) {
  fh <- file.path(chip.external.peak.dir, chip.ext.files[mk])
  peaks.df <- read.csv(fh, header = FALSE, sep = "\t")
  colnames(peaks.df) <- cls2
  peaks.df <- peaks.df[peaks.df[,1] %in% chroms,]
  peaks.gr <- makeGRangesFromDataFrame(peaks.df, keep.extra.columns = TRUE)
  peaks.gr <- sortSeqlevels(peaks.gr)
  peaks.gr <- sort(peaks.gr)
  peaks.gr <- peaks.gr[!overlapsAny(peaks.gr,mm10.blacklist.noovrl.gr)]
  peaks.gr <- annotateRegions(peaks.gr, gencode.txdb, 
                              tss.region = c(-3000,3000),
                              Enh.gr, SEnh.gr, cpg.islands.gr,
                              all_CpG_sites, cpg.meth.gr)
  
  chip.external.peaks.list[[mk]] <- peaks.gr
}


save(list = c("chip.raw.peaks.list", "chip.peaks.list", 
              "chip.raw.rescue.peaks.list",
              "chip.external.peaks.list","chip.rescue.peaks.list",
              "chip.external.stats.df","chip.stats.df",
              "chip.rescue.stats.df"),
     file = file.path(chip.external.data.dir,"ChIP_processed_peaks.RData"))


rm(list = c("chip.raw.peaks.list", "chip.peaks.list", 
            "chip.raw.rescue.peaks.list",
            "chip.external.peaks.list","chip.rescue.peaks.list"))
gc()

#### Compute windows of different sizes of ChIP-seq
#### for various analysis

frag.len <- 150 #round(mean(chip.corr.list[grepl(mark, names(chip.corr.list))]))

spike.chroms <- c("chr2L_spike", "chr2LHet_spike","chr2R_spike", "chr2RHet_spike", 
                  "chr3L_spike", "chr3LHet_spike", "chr3R_spike", "chr3RHet_spike", 
                  "chr4_spike", "chrM_spike", "chrU_spike", "chrUextra_spike",
                  "chrX_spike", "chrXHet_spike","chrYHet_spike")


full.bam.dir <- file.path(chip.data.dir,"MAPPING_mm10_dm3/bam_files")
full.bam.files <- list.files(full.bam.dir,
                             pattern = "bam$")
names(full.bam.files) <- strsplit2(full.bam.files, "\\.")[,1]

input.bam.dir <- file.path(chip.data.dir,"MAPPING_mm10/bam_files")
input.bam.files <-  list.files(input.bam.dir,
                               pattern = ".nodup.bam$")

input.bam.files <- input.bam.files[grepl("input",input.bam.files)]

names(input.bam.files) <- gsub("input", "ChIP",strsplit2(input.bam.files, "\\.")[,1])
input.bam.files <- input.bam.files[names(full.bam.files)]

chip.samples.df <- chip.stats.df
#chip.samples.df <- read.csv(file.path(chip.data.dir,"ChIP_MCM2_TEsubf_sample_info.tsv"),
#                            header = TRUE, sep = "\t")

rownames(chip.samples.df) <- chip.samples.df$ID
chip.samples.df$clone <- factor(chip.samples.df$clone,
                                levels = c("410","421","439","442"))

chip.samples.df$condition <- factor(ifelse(chip.samples.df$clone == "410","WT",
                                           "MCM2_2A"),
                                    levels = c("WT","MCM2_2A"))
chip.samples.df$replicate <- factor(chip.samples.df$replicate)
chip.samples.df$type <- factor(chip.samples.df$type)


window.params.list <- list("w5000" = 
                             list(width = 5000,
                                  spacing = 500),
                           "w1000" = list(width = 1000,
                                          spacing = 100),
                           "w150" = list(width = 150,
                                         spacing = 50))

csaw.cores <- min(16,num.cores)
for (wid in names(window.params.list)) {
  print(wid)
  for (mk in c("H3K27me3","H3K9me3","H3K4me3","H3K27ac","SUZ12")) {
    if (((wid == "w5000") & (mk %in% c("H3K27ac","H3K4me3","SUZ12"))) |
        ((wid == "w150") & (mk %in% c("H3K27me3","H3K9me3","SUZ12")))) {
      next
    }
    print(mk)
    window.width <- window.params.list[[wid]][["width"]]
    window.spacing <- window.params.list[[wid]][["spacing"]]
    
    pe <- "none"
    if (mk %in% (c("H3K4me3", "H3K27ac","SUZ12"))) {
      pe <- "both"
    }
    
    fraglen <- 150
    mids <- c("421","439","442")
    if (mk == "SUZ12") {
      fraglen <- 200
      mids <- c("439","442")
    }
    
    restrict.param.spike <- readParam(discard = mm10.blacklist.noovrl.gr,
                                restrict = c(chroms,spike.chroms), pe = pe)
  
    restrict.param <- readParam(discard = mm10.blacklist.noovrl.gr,
                                      restrict = chroms, pe = pe)
    
    e.bam.files <- full.bam.files[grepl(mk,full.bam.files)]
    i.bam.files <- input.bam.files[names(e.bam.files)]
    
    tmp.samples.df <- filter(chip.samples.df, mark == mk,
                             clone %in% c("410",mids), replicate != "pr")
    
    csaw.chip.bins <- windowCounts(c(file.path(full.bam.dir, e.bam.files), 
                                             file.path(input.bam.dir, i.bam.files)), 
                                           width = 10000, 
                                           bin = TRUE,
                                           ext = frag.len, 
                                           param =  restrict.param, 
                                           BPPARAM = MulticoreParam(workers = csaw.cores))
    
    samp.names <- c(names(e.bam.files), gsub("ChIP","input",names(i.bam.files)))
    
    colnames(csaw.chip.bins) <- samp.names
    
    
    csaw.chip.windows.unspike <- windowCounts(c(file.path(full.bam.dir, e.bam.files), 
                                                file.path(input.bam.dir, i.bam.files)), 
                                            width = window.width, 
                                            spacing = window.spacing,
                                            bin = FALSE,
                                            ext = frag.len,
                                            param = restrict.param,
                                            BPPARAM = MulticoreParam(workers = csaw.cores))
    
    colnames(csaw.chip.windows.unspike) <- samp.names
    
    rownames(csaw.chip.windows.unspike) <- paste0(seqnames(rowRanges(csaw.chip.windows.unspike)),":",
                                               start(rowRanges(csaw.chip.windows.unspike)), "-",
                                               end(rowRanges(csaw.chip.windows.unspike)))

    
    csaw.chip.bins.unspike <- windowCounts(c(file.path(full.bam.dir, e.bam.files), 
                                                file.path(input.bam.dir, i.bam.files)), 
                                              width = window.width, 
                                              spacing = window.spacing,
                                              bin = TRUE,
                                              ext = frag.len,
                                              param = restrict.param,
                                              BPPARAM = MulticoreParam(workers = csaw.cores))
    
    colnames(csaw.chip.bins.unspike) <- samp.names
    
    rownames(csaw.chip.bins.unspike) <- paste0(seqnames(rowRanges(csaw.chip.bins.unspike)),":",
                                                  start(rowRanges(csaw.chip.bins.unspike)), "-",
                                                  end(rowRanges(csaw.chip.bins.unspike)))
    
        
    cids <- c("410",mids)
    csaw.windows.log2FC.mat <- matrix(0,nrow = nrow(csaw.chip.windows.unspike),
                             ncol = length(cids))
    
    rownames(csaw.windows.log2FC.mat) <- rownames(csaw.chip.windows.unspike)
    colnames(csaw.windows.log2FC.mat) <- cids
    
    csaw.bins.log2FC.mat <- matrix(0,nrow = nrow(csaw.chip.bins.unspike),
                                      ncol = length(cids))
    
    rownames(csaw.bins.log2FC.mat) <- rownames(csaw.chip.bins.unspike)
    colnames(csaw.bins.log2FC.mat) <- cids
    
      
    for (cid in cids) {
      chip.samps <- samp.names[grepl(cid, samp.names) & grepl("ChIP",samp.names)]
      input.samps <- samp.names[grepl(cid, samp.names) & grepl("input",samp.names)]
      
      scale.control <- scaleControlFilter(csaw.chip.bins[,chip.samps], 
                                          csaw.chip.bins[,input.samps])
      
      windows.scale.filtered <- filterWindowsControl(csaw.chip.windows.unspike[,chip.samps],
                                             csaw.chip.windows.unspike[,input.samps],
                                             scale.info = scale.control) 
      
      bins.scale.filtered <- filterWindowsControl(csaw.chip.bins.unspike[,chip.samps],
                                                     csaw.chip.bins.unspike[,input.samps],
                                                     scale.info = scale.control) 
      
      
      csaw.windows.log2FC.mat[,cid] <- windows.scale.filtered$filter
      csaw.bins.log2FC.mat[,cid] <- bins.scale.filtered$filter
    }
    
    csaw.chip.windows.spike <- windowCounts(file.path(full.bam.dir, e.bam.files),
                                            width = window.width, 
                                            spacing = window.spacing,
                                            bin = FALSE,
                                            ext = frag.len,
                                            param = restrict.param.spike,
                                            BPPARAM = MulticoreParam(workers = csaw.cores))
    colnames(csaw.chip.windows.spike) <- names(e.bam.files)
    
    rownames(csaw.chip.windows.spike) <- paste0(seqnames(rowRanges(csaw.chip.windows.spike)),":",
                                                  start(rowRanges(csaw.chip.windows.spike)), "-",
                                                  end(rowRanges(csaw.chip.windows.spike)))
    
    
    csaw.chip.windows.exo <- csaw.chip.windows.spike[grepl("_spike",seqnames(rowRanges(csaw.chip.windows.spike))),]
    
    csaw.chip.windows.endo <-  csaw.chip.windows.spike[!grepl("_spike",seqnames(rowRanges(csaw.chip.windows.spike))),]
    
    
    csaw.windows.filter.mat <- csaw.windows.log2FC.mat >= log2(1)
    
    keep.windows <- rownames(csaw.windows.filter.mat)[rowSums(csaw.windows.filter.mat) > 0]
    
    csaw.chip.windows <- csaw.chip.windows.endo[rownames(csaw.chip.windows.endo) %in% keep.windows,]
    
    csaw.windows.log2FC.mat <- csaw.windows.log2FC.mat[rownames(csaw.chip.windows),]
    
    
    tmp.samples.df <- chip.samples.df[colnames(csaw.chip.windows),]
    colData(csaw.chip.windows) <- cbind(colData(csaw.chip.windows),tmp.samples.df)
    out.file <- file.path(chip.data.dir, paste0("ChIP_csaw_windows_",mk ,"_",wid,".RData"))
    save(file = out.file, list = c("csaw.chip.windows",
                                   "csaw.chip.windows.exo",
                                   "csaw.windows.log2FC.mat"))
    
    
    csaw.chip.bins.spike <- windowCounts(file.path(full.bam.dir, e.bam.files),
                                            width = window.width, 
                                            spacing = window.spacing,
                                            bin = FALSE,
                                            ext = frag.len,
                                            param = restrict.param.spike,
                                            BPPARAM = MulticoreParam(workers = csaw.cores))
    colnames(csaw.chip.bins.spike) <- names(e.bam.files)
    
    rownames(csaw.chip.bins.spike) <- paste0(seqnames(rowRanges(csaw.chip.bins.spike)),":",
                                                start(rowRanges(csaw.chip.bins.spike)), "-",
                                                end(rowRanges(csaw.chip.bins.spike)))
    
    
    csaw.chip.bins.exo <- csaw.chip.bins.spike[grepl("_spike",seqnames(rowRanges(csaw.chip.bins.spike))),]
    
    csaw.chip.bins.endo <-  csaw.chip.bins.spike[!grepl("_spike",seqnames(rowRanges(csaw.chip.bins.spike))),]
    
    
    csaw.bins.filter.mat <- csaw.bins.log2FC.mat >= log2(1)
    
    keep.windows <- rownames(csaw.bins.filter.mat)[rowSums(csaw.bins.filter.mat) > 0]
    
    csaw.chip.bins<- csaw.chip.bins.endo[rownames(csaw.chip.bins.endo) %in% keep.windows,]
    
    csaw.bins.log2FC.mat <- csaw.bins.log2FC.mat[rownames(csaw.chip.bins),]
    
    
    tmp.samples.df <- chip.samples.df[colnames(csaw.chip.bins),]
    colData(csaw.chip.bins) <- cbind(colData(csaw.chip.bins),tmp.samples.df)
    out.file <- file.path(chip.data.dir, paste0("ChIP_csaw_bins_",mk ,"_",wid,".RData"))
    save(file = out.file, list = c("csaw.chip.bins",
                                   "csaw.chip.bins.exo",
                                   "csaw.bins.log2FC.mat"))
    
  }
}




    
### Prepare windows for ChIP_rescue analysis

chip.rescue.bam.dir <- file.path(chip.data.dir,
                                 "rescue/MAPPING_mm10_dm3/bam_files")

chip.rescue.bam.files <- list.files(chip.rescue.bam.dir,
                                    pattern = ".bam$")
names(chip.rescue.bam.files) <- strsplit2(chip.rescue.bam.files, "\\.")[,1]

chip.rescue.input.bam.dir <- file.path(chip.data.dir,
                                       "rescue/MAPPING_mm10/bam_files")

chip.rescue.input.bam.files <-  list.files(chip.rescue.input.bam.dir,
                                           pattern = ".nodup.bam$")

names(chip.rescue.input.bam.files) <- gsub("input", "ChIP",strsplit2(chip.rescue.input.bam.files, "\\.")[,1])
chip.rescue.input.bam.files <- chip.rescue.input.bam.files[names(chip.rescue.bam.files)]

csaw.cores <- min(20,num.cores)

window.width <- 1000
window.spacing <- 100
pe <- "none"
frag.len <- 150

restrict.param.spike <- readParam(discard = mm10.blacklist.noovrl.gr,
                                  restrict = c(chroms,spike.chroms), pe = pe)

restrict.param <- readParam(discard = mm10.blacklist.noovrl.gr,
                            restrict = chroms, pe = pe)

e.bam.files <- chip.rescue.bam.files
i.bam.files <- chip.rescue.input.bam.files[names(e.bam.files)]

tmp.samples.df <- chip.rescue.stats.df

csaw.chip.bins <- windowCounts(c(file.path(chip.rescue.bam.dir, e.bam.files), 
                                 file.path(chip.rescue.input.bam.dir, i.bam.files)), 
                               width = 10000, 
                               bin = TRUE,
                               ext = frag.len, 
                               param =  restrict.param, 
                               BPPARAM = MulticoreParam(workers = csaw.cores))

samp.names <- c(names(e.bam.files), gsub("ChIP","input",names(i.bam.files)))

colnames(csaw.chip.bins) <- samp.names


csaw.chip.windows.unspike <- windowCounts(c(file.path(chip.rescue.bam.dir, e.bam.files), 
                                            file.path(chip.rescue.input.bam.dir, i.bam.files)), 
                                          width = window.width, 
                                          spacing = window.spacing,
                                          bin = FALSE,
                                          ext = frag.len,
                                          param = restrict.param,
                                          BPPARAM = MulticoreParam(workers = csaw.cores))

colnames(csaw.chip.windows.unspike) <- samp.names

rownames(csaw.chip.windows.unspike) <- paste0(seqnames(rowRanges(csaw.chip.windows.unspike)),":",
                                              start(rowRanges(csaw.chip.windows.unspike)), "-",
                                              end(rowRanges(csaw.chip.windows.unspike)))

mids <- c("439","442","586","588")
cids <- c("410",mids)
csaw.log2FC.mat <- matrix(0,nrow = nrow(csaw.chip.windows.unspike),
                          ncol = length(cids))

rownames(csaw.log2FC.mat) <- rownames(csaw.chip.windows.unspike)
colnames(csaw.log2FC.mat) <- cids

for (cid in cids) {
  chip.samps <- samp.names[grepl(cid, samp.names) & grepl("ChIP",samp.names)]
  input.samps <- samp.names[grepl(cid, samp.names) & grepl("input",samp.names)]
  
  scale.control <- scaleControlFilter(csaw.chip.bins[,chip.samps], 
                                      csaw.chip.bins[,input.samps])
  
  scale.filtered <- filterWindowsControl(csaw.chip.windows.unspike[,chip.samps],
                                         csaw.chip.windows.unspike[,input.samps],
                                         scale.info = scale.control) 
  csaw.log2FC.mat[,cid] <- scale.filtered$filter
}

csaw.chip.windows.spike <- windowCounts(file.path(chip.rescue.bam.dir, e.bam.files),
                                        width = window.width, 
                                        spacing = window.spacing,
                                        bin = FALSE,
                                        ext = frag.len,
                                        param = restrict.param.spike,
                                        BPPARAM = MulticoreParam(workers = csaw.cores))
colnames(csaw.chip.windows.spike) <- names(e.bam.files)

rownames(csaw.chip.windows.spike) <- paste0(seqnames(rowRanges(csaw.chip.windows.spike)),":",
                                            start(rowRanges(csaw.chip.windows.spike)), "-",
                                            end(rowRanges(csaw.chip.windows.spike)))


csaw.chip.windows.exo <- csaw.chip.windows.spike[grepl("_spike",seqnames(rowRanges(csaw.chip.windows.spike))),]

csaw.chip.windows.endo <-  csaw.chip.windows.spike[!grepl("_spike",seqnames(rowRanges(csaw.chip.windows.spike))),]


csaw.filter.mat <- csaw.log2FC.mat >= log2(1)

keep.windows <- rownames(csaw.filter.mat)[rowSums(csaw.filter.mat) > 0]

csaw.chip.windows <- csaw.chip.windows.endo[rownames(csaw.chip.windows.endo) %in% keep.windows,]

csaw.log2FC.mat <- csaw.log2FC.mat[rownames(csaw.chip.windows),]

out.file <- file.path(chip.data.dir, paste0("ChIP_rescue_csaw_windows_H3K27me3_w1000.RData"))
save(file = out.file, list = c("csaw.chip.windows",
                               "csaw.chip.windows.exo",
                               "csaw.log2FC.mat"))