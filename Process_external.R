library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(org.Mm.eg.db)
library(clusterProfiler)

### Script to process all external data used and in most analyses

wd.dir <- "/isdata/alab/people/nalcaraz/Projects/KU/rep_chromatin_TC/publication/Wenger_et_al_2022"
setwd(wd.dir)

data.dir <- "data"
external.data.dir <- file.path(data.dir,"external")

### Gencode annotations

annot.version <- "gencode.vM23"
assembly <- "mm10"


system("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz")
system(paste("mv gencode.vM23.annotation.gtf.gz",external.data.dir))

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


chroms <- c(paste("chr",1:19,sep=""), "chrX", "chrY")

## Load mm10 blacklist regions and add regiones with deletions
mm10.blacklist.df <- read.csv(file.path(external.data.dir, "mm10.blacklist.bed"),
                              header = FALSE, sep = "\t", 
                              stringsAsFactors = FALSE)
colnames(mm10.blacklist.df) <- c("seqnames", "start", "end")

mm10.blacklist.df$strand <- "*"
mm10.blacklist.gr <- makeGRangesFromDataFrame(mm10.blacklist.df, keep.extra.columns = TRUE)

del.gr <- makeGRangesFromDataFrame(data.frame(seqnames = c("chr6", "chr14"),
                                              start = c(88892208, 47063608),
                                              end = c(91738889, 68368484),
                                              strand = "*", stringsAsFactors = FALSE), keep.extra.columns = TRUE) 

mm10.blacklist.gr <- c(mm10.blacklist.gr, del.gr)
mm10.blacklist.noovrl.gr <- GenomicRanges::reduce(mm10.blacklist.gr, min.gapwidth = 0)

## Load mm10 CpG islands
cpg.islands <- read.csv(file = file.path(external.data.dir,"mm10.cpgIslands.tsv"), sep = "\t", 
                        stringsAsFactors = FALSE, header = FALSE)[,c(2,3,4)]
colnames(cpg.islands) <- c("seqnames", "start", "end")
cpg.islands.gr <- makeGRangesFromDataFrame(cpg.islands, keep.extra.columns = TRUE)

## Load mouse enhancer data in
load(file.path(external.data.dir,"Fantom5_enhancers_robust.RData"))
start(Enh.gr) <- Enh.gr$peakStart
end(Enh.gr) <- Enh.gr$peakEnd
Enh.gr <- Enh.gr[!is.na(Enh.gr$active)]

## Load mouse super enhancers and do liftover
SEnh_meta <- read.csv(file.path(file.path(external.data.dir,"mm9.super_enhancers_meta.bed")),
                      skip = 3,sep=";",header = T)

SEnh <- read.table(file.path(file.path(external.data.dir,
                                       "mm9.super_enhancers.bed")),
                   skip=2,sep=";",header=T)

SEnh.gr <- makeGRangesFromDataFrame(SEnh[SEnh$isSuper=="YES",],seqnames.field = "chrom",
                                    start.field = "start",
                                    end.field = "end",
                                    keep.extra.columns = T)

chain <- import.chain(file.path(external.data.dir,"mm9ToMm10.over.chain"))
SEnh.gr <- unlist(liftOver(SEnh.gr,chain))

ensem2 <- bitr(SEnh.gr$proximal_gene, 
               fromType = "ACCNUM",
               toType = c("ENTREZID"),
               OrgDb = "org.Mm.eg.db")

SEnh.gr$ENTREZID <- ensem2$ENTREZID[match(SEnh.gr$proximal_gene,ensem2$ACCNUM)]
SEnh.gr$name <- paste(seqnames(SEnh.gr),start(SEnh.gr),sep="-") %>% paste(.,end(SEnh.gr),sep=":")
SEnh.gr <- subset(SEnh.gr,seqnames %in% chroms)
SEnh.gr$EPU <- SEnh_meta$Enhancer_and_gene_overlap_EPU[match(SEnh.gr$ID,SEnh_meta$Enhancer_ID)]
SEnh.gr$EPtopo <- SEnh_meta$Enhancer_gene_overlap_topoDomain[match(SEnh.gr$ID,SEnh_meta$Enhancer_ID)]

## Load methylation data

cpg.meth.df <- read.csv(file.path(external.data.dir,"WGBS_mESC_serum_LIF_pool_cpg.txt.gz"),
                        header = TRUE, stringsAsFactors = FALSE, sep = "\t")

cpg.meth.df <- cpg.meth.df[,c(1,2,3,7)]
colnames(cpg.meth.df) <- c("seqnames", "start","end","score")
cpg.meth.gr <- makeGRangesFromDataFrame(cpg.meth.df, keep.extra.columns = TRUE)
cpg.meth.gr$is.methylated <- cpg.meth.gr$score >= 0.7

all_CpG_sites <- vmatchPattern("CG", BSgenome.Mmusculus.UCSC.mm10)
all_CpG_sites <- unstrand(all_CpG_sites[strand(all_CpG_sites) == "+"])

### Load Okazaki-seq
OK.file <- file.path(external.data.dir,"Okazaki_mm10_r1_smooth_results_w1000_s30_d30_z1.txt.gz")

cls <- c("seqnames","start","end","F","R","F.cpm","R.cpm","RFD.raw","RFD",
         "RFD.deriv","score","zero.deriv")

OK.df <- read.table(OK.file,
                    sep="\t",header=FALSE,
                    as.is=TRUE, fill = TRUE)

colnames(OK.df) <- cls
OK.df <- dplyr::filter(OK.df, seqnames %in% chroms)
OK.df$names <- paste(OK.df$seqnames,OK.df$start,sep=":") %>% paste(.,OK.df$end,sep="-") # genomic identifier
OK.df$exprs <- OK.df[,"F"] + OK.df[,"R"] # total raw counts
OK.df$CPM <- OK.df$F.cpm + OK.df$R.cpm # total counts pr. million
OK.gr <- makeGRangesFromDataFrame(OK.df, keep.extra.columns = TRUE)
bl.ovl <- findOverlaps(OK.gr, mm10.blacklist.noovrl.gr, minoverlap = 1,
                       select = "all", type = "any")
OK.gr <- OK.gr[-unique(queryHits(bl.ovl))]
OK.gr.tmp <- subset(OK.gr, CPM >= 0.3)

ids <- which(OK.gr.tmp$zero.deriv > quantile(OK.gr.tmp$RFD.deriv, probs=0.9,
                                             na.rm=TRUE))
OK.gr.tmp <- OK.gr.tmp[ids]
names(ids) <- OK.gr.tmp$names

OK.red.gr <- GenomicRanges::reduce(OK.gr.tmp, min.gapwidth=3000, with.revmap = TRUE)
filtered.data <- OK.gr.tmp[sapply(OK.red.gr$revmap, function(x) {
  x[which.max(OK.gr.tmp$RFD.deriv[x])]
})]

OK.gr$IZ <- ifelse(OK.gr$names %in% filtered.data$names, TRUE, FALSE)
OK.gr$sample <- "OK-seq"

### ======== Load replication timing data ======================================

RT.df <- read.table(file.path(external.data.dir, "RT_binned_results_w1000.txt.gz"),
                    sep="\t",header=FALSE,as.is=TRUE)

RT.df <- RT.df[,c(1:3,6)]
colnames(RT.df) <- c("chr","start","end","log2FC")
RT.df <- filter(RT.df, log2FC != 0)
RT.df$ID <- paste(RT.df$chr,RT.df$start,sep=":") %>% paste(.,RT.df$end,sep="-")
RT.gr <- makeGRangesFromDataFrame(RT.df,
                                  seqnames.field="chr",
                                  start.field = "start",
                                  end.field = "end",
                                  keep.extra.columns = T)
RT.gr <- subset(RT.gr,seqnames %in% chroms)
early.cutoff <- quantile(subset(RT.gr, log2FC > 0)$log2FC, 0.5)
late.cutoff <- quantile(subset(RT.gr, log2FC < 0)$log2FC, 0.5)

RT.gr$RT.class <- factor(sapply(RT.gr$log2FC, FUN = function(x) {
  if (x < 0) {
    if (x < late.cutoff) {
      return("late")  
    }
    return("mid-late")
  } else if (x > early.cutoff) {
    return("early")
  }
  return("mid-early")
  
}), levels = c("early","mid-early","mid-late","late"))

RT.gr <- RT.gr[-unique(queryHits(findOverlaps(RT.gr,mm10.blacklist.noovrl.gr, 
                                              minoverlap = 1, select = "all", type = "any")))]

## Get initiation Zones
IZ.gr <- subset(OK.gr, IZ)

### Save all results

save(file = file.path(external.data.dir, "Processed_external.RData"),
     list = c("all_CpG_sites","chrom.sizes.df","mm10.blacklist.noovrl.gr",
              "cpg.islands.gr","cpg.meth.gr","Enh.gr",
              "SEnh.gr","OK.gr","RT.gr", "IZ.gr"))
