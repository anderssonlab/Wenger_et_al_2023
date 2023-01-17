require(AnnotationDbi)
require(GenomicFeatures)
require(ggalluvial)
require(rtracklayer)
require(CAGEfightR)
require(ggplot2)
theme_set(theme_bw(base_size=6))

## Load txdb object
if (!file.exists("../data/gencode.vM23.annotation.sqlite")) {
    system("wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz")
    system("gzip -d gencode.vM23.annotation.gtf.gz")
    txdb <- makeTxDbFromGFF(file="gencode.vM23.annotation.gtf", format="gtf",
                            chrominfo=read.csv2(file="../data/mm10.chrom.sizes", header=FALSE, col.names=c("chrom", "length"), sep="\t"),
                            organism="Mus musculus")
    saveDb(txdb, file="../data/gencode.vM23.annotation.sqlite")
} else
    txdb <- loadDb("../data/gencode.vM23.annotation.sqlite")

## Retrieve promoters of annotated genes

promoters <- reduce(promoters(txdb,upstream=50,downstream=50))
names(promoters) <- paste0(seqnames(promoters),":",start(promoters),"-",end(promoters),";",strand(promoters))
gene.promoters <- assignTxType(promoters, txdb, tssUpstream = 50, tssDownstream = 50)
gene.promoters <- subset(gene.promoters, txType == "promoter")
gene.promoters <- assignGeneID(gene.promoters, geneModels = txdb, outputColumn = "geneID")
gene.promoters <- gene.promoters[!is.na(gene.promoters$geneID),]
gene.expressed <- scan("../data/genes_background.txt",what=character(0))
gene.promoters.expressed <- gene.promoters[as.character(sapply(mcols(gene.promoters)$geneID,function(n) strsplit(n,"\\.")[[1]][1])) %in% gene.expressed,]

## Chromatin state analysis
H3K4me3 <- read.table("/isdata/alab/projects/rep_chromatin_TC/analysis/data/new_peaks/rescue/410_H3K4me3_ChIP_pr_R1.overlap_peaks.narrowPeak.bed.gz")
H3K27me3 <- read.table("/isdata/alab/projects/rep_chromatin_TC/analysis/data/new_peaks/rescue/410_H3K27me3_ChIP_pr_R1.overlap_50per_peaks.danpos3.regions.bed.gz")
H3K27ac <- read.table("/isdata/alab/projects/rep_chromatin_TC/analysis/data/new_peaks/original/410_H3K27ac_ChIP_pr_R1.overlap_peaks.narrowPeak.bed.gz")
H3K9me3 <- read.table("/isdata/alab/projects/rep_chromatin_TC/analysis/data/new_peaks/rescue/410_H3K9me3_ChIP_pr_R1.overlap_50per_peaks.danpos3.regions.bed.gz")

colnames(H3K4me3) <- colnames(H3K27ac) <- c("chr","start","end","name","score","strand")
colnames(H3K27me3) <- colnames(H3K9me3) <- c("chr","start","end","","strand")

H3K4me3 <- makeGRangesFromDataFrame(H3K4me3)
H3K27me3 <- makeGRangesFromDataFrame(H3K27me3)
H3K27ac <- makeGRangesFromDataFrame(H3K27ac)
H3K9me3 <- makeGRangesFromDataFrame(H3K9me3)

d <- list(H3K4me3=mcols(distanceToNearest(gene.promoters.expressed,H3K4me3,ignore.strand=TRUE))$distance,
          H3K27me3=mcols(distanceToNearest(gene.promoters.expressed,H3K27me3,ignore.strand=TRUE))$distance,
          H3K27ac=mcols(distanceToNearest(gene.promoters.expressed,H3K27ac,ignore.strand=TRUE))$distance,
          H3K9me3=mcols(distanceToNearest(gene.promoters.expressed,H3K9me3,ignore.strand=TRUE))$distance)
gene.promoter.peak.distance <- data.frame(distance=unlist(d),
                                          mark=unlist(lapply(names(d), function(n) rep(n,length(d[[n]])))))

ggplot(gene.promoter.peak.distance,aes(distance+1,color=mark)) + geom_density() + scale_x_log10() + geom_vline(xintercept=c(1000))
ggsave("plots/gene_promoter_nearest_distance_to_marks.pdf",width=8,height=4)


## Create peak set
mcols(H3K4me3)$type <- "H3K4me3"
mcols(H3K27me3)$type <- "H3K27me3"
mcols(H3K27ac)$type <- "H3K27ac"
mcols(H3K9me3)$type <- "H3K9me3"
peaks <- unlist(GRangesList(list(H3K4me3,H3K27me3,H3K27ac,H3K9me3)))

## Assign chromatin states to focus genomic regions according to proximity with ChIP-seq peaks within 1kb

distance <- 1000

gene.promoter.state <- rep("none",length(gene.promoters.expressed))
names(gene.promoter.state) <- names(gene.promoters.expressed)
o <- findOverlaps(gene.promoters.expressed,peaks,maxgap=distance,type="any",select="all")
gene.promoter.state[unique(queryHits(o))] <- by(subjectHits(o),list(queryHits(o)),function(x) paste(sort(unique(peaks$type[x])),collapse=","))
names(gene.promoter.state) <- names(gene.promoters.expressed)

states <- sort(table(as.character(gene.promoter.state)),decreasing=TRUE)

df <- data.frame(promoter=names(gene.promoter.state),
                 gene=as.character(sapply(gene.promoters.expressed[names(gene.promoter.state)]$geneID, function (n) strsplit(n,"\\.")[[1]][1])),
                 state=gene.promoter.state)
write.table(df,file=paste("gene_promoter_states_distance_",distance,".tab",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

## Plot state frequencies

df <- as.data.frame(states)
colnames(df) <- c("State","Frequency")
df$State <- factor(df$State,levels=names(states))
ggplot(df,aes(State,Frequency)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) + theme(axis.text.y = element_text(size=5))
ggsave(paste("plots/Chromatin_state_frequencies_gene_promoters_distance_",distance,".pdf",sep=""),width=6,height=3)

selection <- c(
    "H3K27me3,H3K4me3",
    "H3K27me3,H3K4me3,H3K9me3",
    "H3K27ac,H3K27me3,H3K4me3",
    "H3K27ac,H3K27me3,H3K4me3,H3K9me3",
    "H3K4me3,H3K9me3",
    "H3K27ac,H3K4me3,H3K9me3",
    "H3K27ac,H3K9me3",
    #"H3K27ac,H3K27me3",
    #"H3K27ac,H3K27me3,H3K9me3",
    "H3K27me3",
    "H3K9me3",
    "H3K27me3,H3K9me3",
    "H3K27ac,H3K4me3",
    "H3K27ac",
    "H3K4me3",
    "none")

cat(unique(gene.promoters.expressed[names(which(gene.promoter.state=="H3K27me3,H3K4me3")),]$geneID),file="bivalent_genes.txt",sep="\n")

## Investigate chromatin states with respect to RNA-seq differential expression status

load("../data/RNAseq_DE_genes.RData")

DE.full <- list(down=RNAseq.DE.genes.list[["full"]][["Down"]],
                up=RNAseq.DE.genes.list[["full"]][["Up"]])
DE.442 <- list(down=read.table("/isdata/alab/projects/rep_chromatin_TC/analysis/data/new_peaks/RNAseq_rescue_442_DE_Down.tsv")[-1,1],
               up=read.table("/isdata/alab/projects/rep_chromatin_TC/analysis/data/new_peaks/RNAseq_rescue_442_DE_Up.tsv")[-1,1])
DE.588 <- list(down=read.table("/isdata/alab/projects/rep_chromatin_TC/analysis/data/new_peaks/RNAseq_rescue_588_DE_Down.tsv")[-1,1],
               up=read.table("/isdata/alab/projects/rep_chromatin_TC/analysis/data/new_peaks/RNAseq_rescue_588_DE_Up.tsv")[-1,1])
DE <- list(full=DE.full,
           c442=DE.442,
           c588=DE.588)

theme_set(theme_bw(base_size=9))

gene.promoter.states.DE.enrich <- list()

for (set.name in names(DE)) {
    DE.genes <- DE[[set.name]]
    down <- DE.genes$down
    up <- DE.genes$up
    g.id <- as.character(sapply(gene.promoters.expressed$geneID, function(n) strsplit(n,"\\.")[[1]][1]))
    u <- gene.promoter.state[which(g.id %in% up)]
    d <- gene.promoter.state[which(g.id %in% down)]
    u.g <- g.id[which(g.id %in% up)]
    d.g <- g.id[which(g.id %in% down)]

    df <- rbind(cbind(names(u),u.g,u,"Up"), cbind(names(d),d.g,d,"Down"))
    dimnames(df) <- list(NULL,c("promoter","gene","state","DE"))
    write.table(df,file=paste("gene_promoter_states_DE_",set.name,"_distance_",distance,".tab",sep=""),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

    gene.promoter.states.DE.enrich[[set.name]] <- t(sapply(names(states), function(state) {
        print(state)
        um <- matrix(c(sum(u==state), sum(u!=state), sum(gene.promoter.state==state)-sum(u==state),sum(gene.promoter.state!=state)-sum(u!=state)),byrow=TRUE,ncol=2)
        dm <- matrix(c(sum(d==state), sum(d!=state), sum(gene.promoter.state==state)-sum(d==state),sum(gene.promoter.state!=state)-sum(d!=state)),byrow=TRUE,ncol=2)
        ut <- fisher.test(um)
        dt <- fisher.test(dm)
        c(ut$estimate,ut$p.value,dt$estimate,dt$p.value,sum(u==state),sum(d==state),sum(gene.promoter.state==state))
    }))
    colnames(gene.promoter.states.DE.enrich[[set.name]]) <- c("Up, odds ratio", "Up, p-value", "Down, odds ratio", "Down, p-value","up","down","n")
    rownames(gene.promoter.states.DE.enrich[[set.name]]) <- names(states)

    df <- as.data.frame(rbind(gene.promoter.states.DE.enrich[[set.name]][,c(1,2,5)],
                              gene.promoter.states.DE.enrich[[set.name]][,c(3,4,6)]))
    df$DE <- rep(c("Up","Down"),each=nrow(gene.promoter.states.DE.enrich[[set.name]]))
    colnames(df) <- c("odds.ratio","p.value","count","DE")
    df$state <- rownames(gene.promoter.states.DE.enrich[[set.name]])
    df[df$p.value>0.05,"odds.ratio"] <- NA

    sel <- subset(df, state %in% selection)
    sel$state <- factor(sel$state, levels=rev(selection))

    ggplot(df,aes(x=DE,y=state,fill=log2(odds.ratio))) + geom_tile() +
        scale_fill_distiller(na.value = "grey80", type = "div", palette = "RdBu",
                             rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {scales::rescale_mid(x, to, from, 0)}) +
        geom_text(aes(label=count),size=2) + xlab("") + ylab("") +
        labs(fill="log2(odds ratio)")
    ggsave(paste("plots/gene_promoter_states_distance_",distance,"_enrich_",set.name,"_DE_model_tile.pdf",sep=""),width=4,height=3)

    ggplot(sel,aes(x=DE,y=state,fill=log2(odds.ratio))) + geom_tile() +
        scale_fill_distiller(na.value = "grey80", type = "div", palette = "RdBu",
                             rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {scales::rescale_mid(x, to, from, 0)}) +
        geom_text(aes(label=count),size=2) + xlab("") + ylab("") +
        labs(fill="log2(odds ratio)")
    ggsave(paste("plots/gene_promoter_states_distance_",distance,"_enrich_",set.name,"_DE_model_tile_reorder.pdf",sep=""),width=4,height=3)

    ggplot(sel,aes(x=DE,y=state,fill=log2(odds.ratio))) + geom_tile() +
        scale_fill_distiller(na.value = "grey80", type = "div", palette = "RdBu",
                             rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {scales::rescale_mid(x, to, from, 0)}) +
        labs(fill="log2(odds ratio)")
    ggsave(paste("plots/gene_promoter_states_distance_",distance,"_enrich_",set.name,"_DE_model_tile_reorder_no_numbers.pdf",sep=""),width=4,height=3)
}

## Save session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

