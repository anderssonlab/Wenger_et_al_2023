require(CAGEfightR)
require(AnnotationDbi)
require(GenomicFeatures)
require(ggalluvial)
require(ggplot2)
theme_set(theme_bw(base_size=6))

## Download CAGEfightR extension
system("wget https://github.com/anderssonlab/CAGEfightR_extensions/archive/v0.1.1.tar.gz")
system("tar xvzf v0.1.1.tar.gz")
system("mv CAGEfightR_extensions-0.1.1 CAGEfightR_extensions")

source("CAGEfightR_extensions/decompose.R")


plus_files <- system("ls data/*plus.bw", intern=TRUE)
minus_files <- system("ls data/*minus.bw", intern=TRUE)

names(plus_files) <- names(minus_files) <- as.character(sapply(plus_files, function(n) gsub(".plus.bw","",gsub("CAGE_","",gsub("data/","",n)))))

plus_files <- BigWigFileList(plus_files)
minus_files <- BigWigFileList(minus_files)

## Build CTSS object
orig.CTSSs <- quantifyCTSSs(plusStrand=plus_files,minusStrand=minus_files)
orig.CTSSs <- calcTotalTags(orig.CTSSs)
orig.CTSSs <- calcPooled(orig.CTSSs, inputAssay="counts")
orig.CTSSs <- calcTPM(orig.CTSSs, inputAssay="counts")

## Supported CTSSs are those with at least 2 tags in at least one sample, used for tag clustering
supported.CTSSs <- subsetBySupport(orig.CTSSs,inputAssay="counts",unexpressed=1,minSamples=1)
supported.CTSSs <- calcPooled(supported.CTSSs, inputAssay="counts")

## Identify FANTOM3-stle level 2 tag clusters
orig.TCs <- clusterUnidirectionally(supported.CTSSs)

## Decompose tag clusters by local maxima
decomp.TCs <- decompose(orig.TCs, orig.CTSSs, fn=local_maxima_decompose)

## Quantify tag cluster expression across libraries
decomp.TC.expr <- quantifyClusters(orig.CTSSs, decomp.TCs)
decomp.TC.expr$totalTags <- orig.CTSSs$totalTags
decomp.TC.expr <- calcTPM(decomp.TC.expr,totalTags="totalTags")

orig.TC.expr <- quantifyClusters(orig.CTSSs, orig.TCs)
orig.TC.expr$totalTags <- orig.CTSSs$totalTags
orig.TC.expr <- calcTPM(orig.TC.expr,totalTags="totalTags")

## Build info data frame on libs
condition <- c("410"="WT","439"="MCM2-2A")
info <- data.frame(lib=colnames(decomp.TC.expr),
                   condition=condition[as.character(sapply(colnames(decomp.TC.expr),function(n) strsplit(n,"_")[[1]][1]))],
                   clone=as.character(sapply(colnames(decomp.TC.expr),function(n) strsplit(n,"_")[[1]][1])))
rownames(info) <- info$lib

## Identify TCs supported by expression (>1 tags) in at least two replicates for the same clone
expressed <- assay(decomp.TC.expr,"counts") > 1
clone.expressed <- sapply(unique(info$clone), function(n) rowSums(expressed[,info$clone==n]))
colnames(clone.expressed) <- unique(info$clone)

keep <- rowMax(clone.expressed)>1
supported.TC.expr <- decomp.TC.expr[keep,]
supported.TCs <- decomp.TCs[rownames(supported.TC.expr),]

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

genomeInfo <- bwCommonGenome(plusStrand=plus_files, minusStrand=minus_files, method='intersect')
seqlevels(txdb) <- seqlevels(genomeInfo)

## assign biotype annotations to TCs
supported.TCs <- assignTxType(supported.TCs, txdb, tssUpstream = 500, tssDownstream = 500)
supported.TC.expr <- quantifyClusters(orig.CTSSs, supported.TCs)
supported.TC.expr$totalTags <- orig.CTSSs$totalTags
supported.TC.expr <- calcTPM(supported.TC.expr,totalTags="totalTags")

## Retrieve promoters of annotated genes
promoters <- reduce(promoters(txdb,upstream=50,downstream=50))
names(promoters) <- paste0(seqnames(promoters),":",start(promoters),"-",end(promoters),";",strand(promoters))
gene.promoters <- assignTxType(promoters, txdb, tssUpstream = 50, tssDownstream = 50)
gene.promoters <- subset(gene.promoters, txType == "promoter")
gene.promoters <- assignGeneID(gene.promoters, geneModels = txdb, outputColumn = "geneID")
gene.promoters <- gene.promoters[!is.na(gene.promoters$geneID),]

## Filter proximal TCs by expression prioritization
tc.olap <- findOverlaps(supported.TCs,supported.TCs,maxgap=250,type="any")
tc.olap <- tc.olap[-which(subjectHits(tc.olap)==queryHits(tc.olap))]

o <- order(supported.TCs$score,decreasing=TRUE)
rem <- c()

for (i in o) {
    if (i %in% rem)
        next
    rem <- c(rem,subjectHits(tc.olap)[which(queryHits(tc.olap)==i)])
}
rem <- unique(rem)

standalone.TCs <- supported.TCs[-rem,]
standalone.TC.expr <- supported.TC.expr[names(standalone.TCs),]


## Filter external focus TSSs by gene expression support
gene.expressed <- scan("../data/genes_background.txt",what=character(0))
gene.promoters.expressed <- gene.promoters[as.character(sapply(mcols(gene.promoters)$geneID,function(n) strsplit(n,"\\.")[[1]][1])) %in% gene.expressed,]


## Chromatin state analysis
load("../data/ChIP_all_peaks.RData")

d <- list(H3K4me3=mcols(distanceToNearest(gene.promoters.expressed,chip.k4.clone.peaks.grl[["410"]],ignore.strand=TRUE))$distance,
          H3K27me3=mcols(distanceToNearest(gene.promoters.expressed,chip.k27.clone.peaks.grl[["410"]],ignore.strand=TRUE))$distance,
          H3K27ac=mcols(distanceToNearest(gene.promoters.expressed,chip.k27ac.clone.peaks.grl[["410"]],ignore.strand=TRUE))$distance,
          H3K9me3=mcols(distanceToNearest(gene.promoters.expressed,chip.k9.clone.peaks.grl[["410"]],ignore.strand=TRUE))$distance)
gene.promoter.peak.distance <- data.frame(distance=unlist(d),
                                          mark=unlist(lapply(names(d), function(n) rep(n,length(d[[n]])))))

ggplot(gene.promoter.peak.distance,aes(distance+1,color=mark)) + geom_density() + scale_x_log10() + geom_vline(xintercept=c(1000))
ggsave("plots/gene_promoter_nearest_distance_to_marks.pdf",width=8,height=4)


## Create peak set

clones <- c("410","421","439","442")

peaks <- lapply(clones, function(clone) {
    print(clone)

    K4 <- chip.k4.clone.peaks.grl[[clone]]
    mcols(K4)$type <- "H3K4me3"

    K27 <- chip.k27.clone.peaks.grl[[clone]]
    mcols(K27)$type <- "H3K27me3"

    K9 <- chip.k9.clone.peaks.grl[[clone]]
    mcols(K9)$type <- "H3K9me3"

    K27ac <- chip.k27ac.clone.peaks.grl[[clone]]
    mcols(K27ac)$type <- "H3K27ac"

    unlist(GRangesList(list(K4,K27,K9,K27ac)))
})
names(peaks) <- clones
peaks <- GRangesList(peaks)


## Assign chromatin states to focus genomic regions according to proximity with ChIP-seq peaks within 1kb

distance <- 1000

standalone.TC.states <- lapply(clones, function(clone) {
    print(clone)
    state <- rep("none",length(standalone.TCs))
    names(state) <- names(standalone.TCs)
    o <- findOverlaps(standalone.TCs,peaks[[clone]],maxgap=distance,type="any",select="all")
    t <- peaks[[clone]]$type
    state[unique(queryHits(o))] <- by(subjectHits(o),list(queryHits(o)),function(x) paste(sort(unique(t[x])),collapse=","))
    state
})
names(standalone.TC.states) <- clones

gene.promoter.states <- lapply(clones, function(clone) {
    print(clone)
    state <- rep("none",length(gene.promoters.expressed))
    names(state) <- names(gene.promoters.expressed)
    o <- findOverlaps(gene.promoters.expressed,peaks[[clone]],maxgap=distance,type="any",select="all")
    t <- peaks[[clone]]$type
    state[unique(queryHits(o))] <- by(subjectHits(o),list(queryHits(o)),function(x) paste(sort(unique(t[x])),collapse=","))
    state
})
names(gene.promoter.states) <- clones

states <- sort(table(c(as.character(unlist(standalone.TC.states)),as.character(unlist(gene.promoter.states)))),decreasing=TRUE)

## Plot state frequencies

df <- do.call("rbind",lapply(names(standalone.TC.states),function(clone) {
                          x <- cbind(as.data.frame(table(standalone.TC.states[[clone]])),clone)
                          colnames(x) <- c("State","Freq","Clone")
                          x
                      }))
df$State <- factor(df$State,levels=names(states))
ggplot(df,aes(State,Clone,fill=Freq)) + geom_tile() + scale_fill_distiller(type = "seq", palette = "Oranges", direction=1, trans = 'log10') +
    geom_text(aes(label=Freq),size=1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) + theme(axis.text.y = element_text(size=5))
ggsave("plots/Chromatin_state_frequencies_standalone_TCs.pdf",width=10,height=5)

df <- do.call("rbind",lapply(names(gene.promoter.states),function(clone) {
                          x <- cbind(as.data.frame(table(gene.promoter.states[[clone]])),clone)
                          colnames(x) <- c("State","Freq","Clone")
                          x
                      }))
df$State <- factor(df$State,levels=names(states))
ggplot(df,aes(State,Clone,fill=Freq)) + geom_tile() + scale_fill_distiller(type = "seq", palette = "Oranges", direction=1, trans = 'log10') +
    geom_text(aes(label=Freq),size=1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) + theme(axis.text.y = element_text(size=5))
ggsave("plots/Chromatin_state_frequencies_gene_promoters.pdf",width=10,height=5)


standalone.TC.states.transitions <- do.call("rbind",
                                            lapply(c("421","439","442"),function(n) {
                                                freq <- table(as.data.frame(cbind(standalone.TC.states[["410"]],
                                                                                  standalone.TC.states[[n]])))
                                                m <- matrix(0,nrow=length(states),ncol=length(states),
                                                            dimnames=list(names(states),names(states)))
                                                m[rownames(freq),colnames(freq)] <- as.matrix(freq)
                                                x  <- data.frame(WT=rownames(m),
                                                                 MUT=rep(colnames(m),each=nrow(m)),
                                                                 Freq=as.numeric(m),
                                                                 sample=n)
                                                x$fraction <- (x$Freq+1) / sum(x$Freq+1)
                                                x
                                            }))

gene.promoter.states.transitions <- do.call("rbind",
                                            lapply(c("421","439","442"),function(n) {
                                                freq <- table(as.data.frame(cbind(gene.promoter.states[["410"]],
                                                                                  gene.promoter.states[[n]])))
                                                m <- matrix(0,nrow=length(states),ncol=length(states),
                                                            dimnames=list(names(states),names(states)))
                                                m[rownames(freq),colnames(freq)] <- as.matrix(freq)
                                                x  <- data.frame(WT=rownames(m),
                                                                 MUT=rep(colnames(m),each=nrow(m)),
                                                                 Freq=as.numeric(m),
                                                                 sample=n)
                                                x$fraction <- (x$Freq+1) / sum(x$Freq+1)
                                                x
                                            }))


## Order states by similarities
df.mean <- subset(standalone.TC.states.transitions,sample=="439")[,c("WT","MUT","Freq","fraction")]
mat.mean <- matrix(df.mean$fraction,nrow=length(unique(df.mean$WT)),dimnames=list(unique(df.mean$WT),unique(df.mean$WT)),byrow=TRUE)
o <- hclust(dist(mat.mean))$order
standalone.TC.states.transitions$WT <- factor(standalone.TC.states.transitions$WT, levels=rownames(mat.mean)[o])
standalone.TC.states.transitions$MUT <- factor(standalone.TC.states.transitions$MUT, levels=rownames(mat.mean)[o])

df.mean <- subset(gene.promoter.states.transitions,sample=="439")[,c("WT","MUT","Freq","fraction")]
mat.mean <- matrix(df.mean$fraction,nrow=length(unique(df.mean$WT)),dimnames=list(unique(df.mean$WT),unique(df.mean$WT)),byrow=TRUE)
o <- hclust(dist(mat.mean))$order
gene.promoter.states.transitions$WT <- factor(gene.promoter.states.transitions$WT, levels=rownames(mat.mean)[o])
gene.promoter.states.transitions$MUT <- factor(gene.promoter.states.transitions$MUT, levels=rownames(mat.mean)[o])

## Plot transitions

ggplot(standalone.TC.states.transitions,aes(x=MUT,y=WT,fill=fraction)) + geom_tile() + scale_fill_distiller(type = "seq", palette = "Oranges", direction=1, trans = 'log10') +
    geom_text(aes(label=Freq),size=2) + xlab("MUT") + ylab("WT") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(sample ~ ., ncol=3)
ggsave("plots/Chromatin_state_tile_standalone_TC_states.pdf",width=18,height=5)

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

sel <- subset(standalone.TC.states.transitions, WT %in% selection & MUT %in% selection)
sel$WT <- factor(sel$WT, levels=rev(selection))
sel$MUT <- factor(sel$MUT, levels=rev(selection))

ggplot(sel,aes(x=MUT,y=WT,fill=fraction)) + geom_tile() + scale_fill_distiller(type = "seq", palette = "Oranges", direction=1, trans = 'log10') +
geom_text(aes(label=Freq),size=2) + xlab("MUT") + ylab("WT") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(sample ~ ., ncol=3)
ggsave("plots/Chromatin_state_tile_standalone_TC_states_reorder.pdf",width=18,height=5)

ggplot(gene.promoter.states.transitions,aes(x=MUT,y=WT,fill=fraction)) + geom_tile() + scale_fill_distiller(type = "seq", palette = "Oranges", direction=1, trans = 'log10') +
    geom_text(aes(label=Freq),size=2) + xlab("MUT") + ylab("WT") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(sample ~ ., ncol=3)
ggsave("plots/Chromatin_state_tile_gene_promoter_states.pdf",width=18,height=5)

ggplot(standalone.TC.states.transitions, aes(y = Freq, axis1 = WT, axis2 = MUT)) +
    geom_alluvium(aes(fill = WT)) +
    geom_stratum(fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("WT (410)", "MUT"), expand = c(.05, .05)) +
    facet_wrap(sample ~ ., ncol=3) +
    ggtitle("Chromatin state")
ggsave("plots/Chromatin_state_alluvial_standalone_TC_states.pdf",width=32,height=16)

ggplot(gene.promoter.states.transitions, aes(y = Freq, axis1 = WT, axis2 = MUT)) +
    geom_alluvium(aes(fill = WT)) +
    geom_stratum(fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("WT (410)", "MUT"), expand = c(.05, .05)) +
    facet_wrap(sample ~ ., ncol=3) +
    ggtitle("Chromatin state")
ggsave("plots/Chromatin_state_alluvial_gene_promoter_states.pdf",width=32,height=16)

## Print bivalent genes for each clone
for (clone in clones)
    cat(unique(gene.promoters[names(which(gene.promoter.states[[clone]]=="H3K27me3,H3K4me3")),]$geneID),file=paste("bivalent_genes_",clone,".txt",sep=""),sep="\n")

## Investigate chromatin state transitions with respect to RNA-seq differential expression status

load("../data/RNAseq_DE_genes.RData")
down <- RNAseq.DE.genes.list[["full"]][["Down"]]
up <- RNAseq.DE.genes.list[["full"]][["Up"]]

wt <- gene.promoter.states[["410"]]
g.id <- as.character(sapply(gene.promoters.expressed$geneID, function(n) strsplit(n,"\\.")[[1]][1]))
u <- wt[which(g.id %in% up)]
d <- wt[which(g.id %in% down)]

u.g <- g.id[which(g.id %in% up)]
d.g <- g.id[which(g.id %in% down)]

df <- rbind(cbind(names(u),u.g,u,"Up"), cbind(names(d),d.g,d,"Down"))
dimnames(df) <- list(NULL,c("promoter","gene","state","DE"))
write.table(df,file="gene_promoter_states_WT.tab",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

gene.promoter.states.WT.DE.full.enrich <- t(sapply(names(states), function(state) {
    print(state)
    um <- matrix(c(sum(u==state), sum(u!=state), sum(wt==state)-sum(u==state),sum(wt!=state)-sum(u!=state)),byrow=TRUE,ncol=2)
    dm <- matrix(c(sum(d==state), sum(d!=state), sum(wt==state)-sum(d==state),sum(wt!=state)-sum(d!=state)),byrow=TRUE,ncol=2)
    ut <- fisher.test(um)
    dt <- fisher.test(dm)
    c(ut$estimate,ut$p.value,dt$estimate,dt$p.value,sum(u==state),sum(d==state),sum(wt==state))
}))
colnames(gene.promoter.states.WT.DE.full.enrich) <- c("Up, odds ratio", "Up, p-value", "Down, odds ratio", "Down, p-value","up","down","n")
rownames(gene.promoter.states.WT.DE.full.enrich) <- names(states)


theme_set(theme_bw(base_size=9))

df <- as.data.frame(rbind(gene.promoter.states.WT.DE.full.enrich[,c(1,2,5)],
                          gene.promoter.states.WT.DE.full.enrich[,c(3,4,6)]))
df$DE <- rep(c("Up","Down"),each=nrow(gene.promoter.states.WT.DE.full.enrich))
colnames(df) <- c("odds.ratio","p.value","count","DE")
df$state <- rownames(gene.promoter.states.WT.DE.full.enrich)
df[df$p.value>0.05,"odds.ratio"] <- NA

sel <- subset(df, state %in% selection)
sel$state <- factor(sel$state, levels=rev(selection))

ggplot(df,aes(x=DE,y=state,fill=log2(odds.ratio))) + geom_tile() +
    scale_fill_distiller(na.value = "grey80", type = "div", palette = "RdBu",
                         rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {scales::rescale_mid(x, to, from, 0)}) +
    geom_text(aes(label=count),size=2) + xlab("") + ylab("")
ggsave("plots/gene_promoter_states_WT_enrich_full_DE_model_tile.pdf",width=4,height=4)

ggplot(sel,aes(x=DE,y=state,fill=log2(odds.ratio))) + geom_tile() +
    scale_fill_distiller(na.value = "grey80", type = "div", palette = "RdBu",
                         rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {scales::rescale_mid(x, to, from, 0)}) +
    geom_text(aes(label=count),size=2) + xlab("") + ylab("")
ggsave("plots/gene_promoter_states_WT_enrich_full_DE_model_tile_reorder.pdf",width=4,height=4)

ggplot(sel,aes(x=DE,y=state,fill=log2(odds.ratio))) + geom_tile() +
    scale_fill_distiller(na.value = "grey80", type = "div", palette = "RdBu",
                         rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {scales::rescale_mid(x, to, from, 0)})
ggsave("plots/gene_promoter_states_WT_enrich_full_DE_model_tile_reorder_no_numbers.pdf",width=4,height=4)


g.id <- as.character(sapply(gene.promoters.expressed$geneID, function(n) strsplit(n,"\\.")[[1]][1]))
u.id <- which(g.id %in% up)
d.id <- which(g.id %in% down)
gene.promoter.states.transitions.DE <- do.call("rbind",
                                               lapply(c("421","439","442"),function(n) {
                                                   ufreq <- table(as.data.frame(cbind(gene.promoter.states[["410"]][u.id],
                                                                                      gene.promoter.states[[n]][u.id])))
                                                   dfreq <- table(as.data.frame(cbind(gene.promoter.states[["410"]][d.id],
                                                                                      gene.promoter.states[[n]][d.id])))
                                                   um <- matrix(0,nrow=length(states),ncol=length(states),
                                                                dimnames=list(names(states),names(states)))
                                                   dm <- um
                                                   um[rownames(ufreq),colnames(ufreq)] <- as.matrix(ufreq)
                                                   dm[rownames(dfreq),colnames(dfreq)] <- as.matrix(dfreq)
                                                   xu  <- data.frame(WT=rownames(um),
                                                                     MUT=rep(colnames(um),each=nrow(um)),
                                                                     Freq=as.numeric(um),
                                                                     sample=n,
                                                                     DE="Up")
                                                   xd  <- data.frame(WT=rownames(dm),
                                                                     MUT=rep(colnames(dm),each=nrow(dm)),
                                                                     Freq=as.numeric(dm),
                                                                     sample=n,
                                                                     DE="Down")
                                                   xu$fraction <- (xu$Freq+1) / sum(xu$Freq+1)
                                                   xd$fraction <- (xd$Freq+1) / sum(xd$Freq+1)
                                                   rbind(xu,xd)
                                               }))


gene.promoter.states.transitions.DE.full.enrich <- do.call("rbind",
                                                           lapply(c("421","439","442"),function(n) {
                                                               print(n)
                                                               x <- subset(gene.promoter.states.transitions.DE,sample==n)
                                                               u <- subset(x,DE=="Up")
                                                               d <- subset(x,DE=="Down")
                                                               a <- subset(gene.promoter.states.transitions,sample==n)
                                                               res <-
                                                                   do.call("rbind",
                                                                           lapply(names(states), function(state1) {
                                                                               t(sapply(names(states), function(state2) {
                                                                                   um <- matrix(c(subset(u,WT==state1 & MUT==state2)$Freq,
                                                                                                  sum(u$Freq)-subset(u,WT==state1 & MUT==state2)$Freq,
                                                                                                  subset(a,WT==state1 & MUT==state2)$Freq-subset(u,WT==state1 & MUT==state2)$Freq,
                                                                                   (sum(a$Freq)-subset(a,WT==state1 & MUT==state2)$Freq)-
                                                                                   (sum(u$Freq)-subset(u,WT==state1 & MUT==state2)$Freq)),
                                                                                   byrow=TRUE,ncol=2)
                                                                                   dm <- matrix(c(subset(d,WT==state1 & MUT==state2)$Freq,
                                                                                                  sum(d$Freq)-subset(d,WT==state1 & MUT==state2)$Freq,
                                                                                                  subset(a,WT==state1 & MUT==state2)$Freq-subset(d,WT==state1 & MUT==state2)$Freq,
                                                                                   (sum(a$Freq)-subset(a,WT==state1 & MUT==state2)$Freq)-
                                                                                   (sum(d$Freq)-subset(d,WT==state1 & MUT==state2)$Freq)),
                                                                                   byrow=TRUE,ncol=2)
                                                                                   ut <- fisher.test(um)
                                                                                   dt <- fisher.test(dm)
                                                                                   c(state1, state2, n,
                                                                                     ut$estimate,ut$p.value,dt$estimate,dt$p.value,
                                                                                     subset(u,WT==state1 & MUT==state2)$Freq,
                                                                                     subset(d,WT==state1 & MUT==state2)$Freq,
                                                                                     subset(a,WT==state1 & MUT==state2)$Freq)
                                                                               }))
                                                                           }))
                                                               res <- as.data.frame(res)
                                                               for (i in 4:10)
                                                                   res[,i] <- as.numeric(res[,i])
                                                               res
                                                           }))
colnames(gene.promoter.states.transitions.DE.full.enrich) <- c("WT","MUT","sample","up.odds.ratio", "up.p.value", "down.odds.ratio", "down.p.value","up","down","n")


df <- gene.promoter.states.transitions.DE.full.enrich
df[df$up.p.value>0.05,"up.odds.ratio"] <- NA
df[df$down.p.value>0.05,"down.odds.ratio"] <- NA

ggplot(df,aes(x=MUT,y=WT,fill=log2(up.odds.ratio))) + geom_tile() +
    scale_fill_distiller(type = "div", palette = "RdBu",
                         rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {scales::rescale_mid(x, to, from, 0)}) +
    geom_text(aes(label=up),size=2) + xlab("") + ylab("") + facet_wrap(~sample,ncol=3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/gene_promoter_states_transitions_enrich_full_DE_Up_model_tile.pdf",width=10,height=4)

ggplot(df,aes(x=MUT,y=WT,fill=log2(down.odds.ratio))) + geom_tile() +
    scale_fill_distiller(type = "div", palette = "RdBu",
                         rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {scales::rescale_mid(x, to, from, 0)}) +
    geom_text(aes(label=down),size=2) + xlab("") + ylab("") + facet_wrap(~sample,ncol=3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/gene_promoter_states_transitions_enrich_full_DE_Down_model_tile.pdf",width=10,height=4)


## Save session info
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
