### Utilty functions for Wenger et al. R scripts 
### Nicolas Alcaraz (2019)


### Function to annotate genomic regions with 
### Extension from ChIPseeker's annotatePeak function
### can give enhancer and cpg.island regions

annotateRegions <- function(granges.obj, txdb.obj, 
                            tss.region = c(-3000,3000),
                            enhancers, superenhancers, cpg.islands,
                            all_CpG_sites, cpg.meth.gr) {
  require(ChIPseeker)
  require(org.Mm.eg.db)
  require(BSgenome.Mmusculus.UCSC.mm10)
  require(GenomicRanges)
  require(GenomicFeatures)
  require(Biostrings)  
  
  names(granges.obj) <- paste(seqnames(granges.obj), ":", start(granges.obj), "-", end(granges.obj), sep = "")
  granges.obj$name <- names(granges.obj)  
  annot <- annotatePeak(granges.obj, tssRegion = tss.region, 
                             TxDb = txdb.obj, annoDb = "org.Mm.eg.db",
                             level = "transcript", sameStrand = FALSE)
  
  # currate UCSC annotations
  annot.df <- data.frame(annot)
  annot.df$annotation_cur <- annot.df$annotation
  annot.df$annotation_cur[grep("exon",annot.df$annotation)] <- "exon"
  annot.df$annotation_cur[grep("intron",annot.df$annotation)] <- "intron"
  annot.df$annotation_cur[grep("Downstream",annot.df$annotation)] <- "Downstream <3kb"
  annot.df$annotation_cur[grep("Promoter",annot.df$annotation)] <- "Promoter"
  
  
  # add annotations to data range
  granges.obj$annotation <- annot.df$annotation
  granges.obj$annotation_cur <- annot.df$annotation_cur
  granges.obj$SYMBOL <- annot.df$SYMBOL
  granges.obj$GENENAME <- annot.df$GENENAME
  granges.obj$ENTREZID <- annot.df$ENTREZID
  granges.obj$ENSEMBL <- annot.df$geneId
  granges.obj$ENSEMBL.transcript <- annot.df$transcriptId
  
  granges.obj$geneStrand <- annot.df$geneStrand
  granges.obj$geneStrand <- ifelse(annot.df$geneStrand=="1","+","-")
  granges.obj$distToTSS <- annot.df$distanceToTSS
  granges.obj$geneLength <- annot.df$geneLength
  granges.obj$geneStart <- annot.df$geneStart
  granges.obj$geneEnd <- annot.df$geneEnd
 
  enh_ovarlap <- findOverlaps(granges.obj, enhancers)
  
  granges.obj$enhancer <- FALSE
  granges.obj$enhancerNo <- FALSE
  if (!missing(enhancers)) {
      granges.obj$enhancer[queryHits(enh_ovarlap)] <-  Enh.gr$active[subjectHits(enh_ovarlap)]
      granges.obj$enhancerNo[queryHits(enh_ovarlap)] <-  Enh.gr$name[subjectHits(enh_ovarlap)]
  }
      
  granges.obj$Senhancer <- FALSE  
  
  if (!missing(superenhancers)) {
      Senh_ovarlap <- findOverlaps(granges.obj, superenhancers)
      granges.obj$Senhancer[queryHits(Senh_ovarlap)] <-  TRUE
  }
    
  granges.obj$enhancer.class <- "None"
  if (sum(granges.obj$enhancer %in% c("active","inter")) > 0) {
    granges.obj[granges.obj$enhancer %in% c("active","inter")]$enhancer.class <- "Active enhancer"
  }
  if (sum(granges.obj$enhancer == "inactive") > 0) {
    granges.obj[granges.obj$enhancer == "inactive"]$enhancer.class <- "Inactive enhancer"
  }
  if (sum(granges.obj$Senhancer) > 0) {
    granges.obj[granges.obj$Senhancer]$enhancer.class <- "Super enhancer"
  }
  granges.obj$enhancer.class <- factor(granges.obj$enhancer.class, 
                                                  levels = c("Super enhancer", "Active enhancer", "Inactive enhancer", "None"))
  
  granges.obj.freq <- alphabetFrequency(getSeq(BSgenome.Mmusculus.UCSC.mm10, granges.obj), baseOnly = TRUE)
  granges.obj$GC.cont <- rowSums(granges.obj.freq[,c("G","C")]) / rowSums(granges.obj.freq[,c("G","C", "A", "T")])
  granges.obj$CpG.island <- overlapsAny(granges.obj, cpg.islands)
  granges.obj$CpG.number <- countOverlaps(granges.obj, all_CpG_sites, minoverlap = 1)
  granges.obj$CpG.content <- 0
  granges.obj$CpG.content <- 0
  granges.obj$GC.cont[is.na(granges.obj$GC.cont)] <- 0.5
  indx <- granges.obj$GC.cont > 0 #& !is.na(granges.obj$GC.cont)
  granges.obj$CpG.content[indx] <- 
    granges.obj$CpG.number[indx] / (((granges.obj.freq[indx,"G"]) * (granges.obj.freq[indx,"C"])) / width(granges.obj[indx]))
  
  granges.obj$mCpG.number <- countOverlaps(granges.obj, subset(cpg.meth.gr, is.methylated), 
                                           minoverlap = 1)
  
  granges.obj$mCpG.fraction <- 0
  granges.obj$mCpG.fraction[granges.obj$CpG.number > 0] <-  granges.obj[granges.obj$CpG.number > 0]$mCpG.number / granges.obj[granges.obj$CpG.number > 0]$CpG.number
  
  granges.obj$mCpG.content <- 0
  granges.obj$mCpG.content <- 0
  indx <- granges.obj$GC.cont > 0
  granges.obj$mCpG.content[indx] <- 
    granges.obj$mCpG.number[indx] / (((granges.obj.freq[indx,"G"]) * (granges.obj.freq[indx,"C"])) / width(granges.obj[indx]))
  
  granges.obj$CpG.content[is.na(granges.obj$CpG.content)] <- 0
  granges.obj$mCpG.content[is.na(granges.obj$mCpG.content)] <- 0
  
  dist.iz <- distanceToNearest(resize(granges.obj, fix = "center", width = 1),
                               resize(IZ.gr, fix = "center", width = 1))
  
  granges.obj$distToIZ <- start(resize(granges.obj[queryHits(dist.iz)], fix = "center", width = 1)) -
    start(resize(IZ.gr[subjectHits(dist.iz)], fix = "center", width = 1))
  granges.obj$distToIZnorm <- 1
  for (chr in chroms) {
    tmp.dist <- distanceToNearest(resize(IZ.gr[seqnames(IZ.gr) == chr], fix = "center", width = 1))
    mx.dist <- max(as.data.frame(tmp.dist)$distance) / 2
    granges.obj[seqnames(granges.obj) == chr]$distToIZnorm <-
      granges.obj[seqnames(granges.obj) == chr]$distToIZ / mx.dist
  }
  granges.obj$distToIZnorm[granges.obj$distToIZnorm > 1] <- 1
  granges.obj$distToIZnorm[granges.obj$distToIZnorm < -1] <- -1
  
 
  return(granges.obj)
}

### Figures for counts
pval2signf <- function(x) {
  return(sapply(x, FUN = function(p) {
    if (p < 0.1 & p >= 0.01) {
      return(".")
    } else if (p < 0.01 & p >= 0.001) {
      return("*")
    } else if (p < 0.001 & p >= 0.0001) {
      return("**") 
    } else if (p < 0.0001 & p >= 0.00001) {
      return("***") 
    } else if (p < 0.00001) {
      return("****")
    }
    return("")  
  }))
}


