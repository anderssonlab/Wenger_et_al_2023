require(CAGEfightR)
require(AnnotationDbi)
require(GenomicFeatures)
require(ggalluvial)
require(ggplot2)
theme_set(theme_bw(base_size=10))

## Download CAGEfightR extension
system("wget https://github.com/anderssonlab/CAGEfightR_extensions/archive/v0.1.0.tar.gz")
system("tar xvzf v0.1.0.tar.gz")
system("mv CAGEfightR_extensions-0.1.0 CAGEfightR_extensions")

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
TCs <- decompose(orig.TCs, orig.CTSSs, fn=local_maxima_decompose)

## Quantify tag cluster expression across libraries
TC.expr <- quantifyClusters(orig.CTSSs, TCs)
TC.expr$totalTags <- orig.CTSSs$totalTags
TC.expr <- calcTPM(TC.expr,totalTags="totalTags")

orig.TC.expr <- quantifyClusters(orig.CTSSs, orig.TCs)
orig.TC.expr$totalTags <- orig.CTSSs$totalTags
orig.TC.expr <- calcTPM(orig.TC.expr,totalTags="totalTags")

## Build info data frame on libs
condition <- c("410"="WT","506"="WT","438"="MCM2-2A","439"="MCM2-2A")
info <- data.frame(lib=colnames(TC.expr),
                   condition=condition[as.character(sapply(colnames(TC.expr),function(n) strsplit(n,"_")[[1]][1]))],
                   clone=as.character(sapply(colnames(TC.expr),function(n) strsplit(n,"_")[[1]][1])))
rownames(info) <- info$lib

## Identify TCs supported by expression (>1 tags) in at least two replicates for the same clone
expressed <- assay(TC.expr,"counts") > 1
clone.expressed <- sapply(unique(info$clone), function(n) rowSums(expressed[,info$clone==n]))
colnames(clone.expressed) <- unique(info$clone)

keep <- rowMax(clone.expressed)>1
supported.TC.expr <- TC.expr[keep,]
supported.TCs <- TCs[rownames(supported.TC.expr),]

## Load txdb object
if (!file.exists("../data/gencode.vM24.annotation.sqlite")) {
    system("wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz")
    system("gzip -d gencode.vM24.annotation.gtf.gz")
    txdb <- makeTxDbFromGFF(file="gencode.vM24.annotation.gtf", format="gtf",
                            chrominfo=read.csv2(file="../data/mm10.chrom.sizes", header=FALSE, col.names=c("chrom", "length"), sep="\t"),
                            organism="Mus musculus")
    saveDb(txdb, file="../data/gencode.vM24.annotation.sqlite")
} else
    txdb <- loadDb("../data/gencode.vM24.annotation.sqlite")

genomeInfo <- bwCommonGenome(plusStrand=plus_files, minusStrand=minus_files, method='intersect')
seqlevels(txdb) <- seqlevels(genomeInfo)

## assign biotype annotations to TCs
supported.TCs <- assignTxType(supported.TCs, txdb, tssUpstream = 500, tssDownstream = 500)
supported.TC.expr <- quantifyClusters(orig.CTSSs, supported.TCs)
supported.TC.expr$totalTags <- orig.CTSSs$totalTags
supported.TC.expr <- calcTPM(supported.TC.expr,totalTags="totalTags")

## Filter TCs proximal to annotated gene TSSs
gene.TCs <- assignTxType(supported.TCs, txdb, tssUpstream = 500, tssDownstream = 500)
gene.TCs <- subset(gene.TCs, txType == "promoter")
gene.TCs <- assignGeneID(gene.TCs, geneModels = txdb, outputColumn = "geneID")
gene.TC.expr <- quantifyClusters(orig.CTSSs, gene.TCs)
gene.TC.expr$totalTags <- orig.CTSSs$totalTags
gene.TC.expr <- calcTPM(gene.TC.expr,totalTags="totalTags")
gene.expr <- quantifyGenes(gene.TC.expr,"geneID")

## Filter proximal TCs
standalone.TCs <- supported.TCs
tc.olap <- findOverlaps(standalone.TCs,standalone.TCs,maxgap=500,type="any")
tc.olap <- tc.olap[-which(subjectHits(tc.olap)==queryHits(tc.olap))]

o <- order(standalone.TCs$score,decreasing=TRUE)
rem <- c()

for (i in o) {
    if (i %in% rem)
        next
    rem <- c(rem,subjectHits(tc.olap)[which(queryHits(tc.olap)==i)])
}
rem <- unique(rem)

standalone.TCs <- standalone.TCs[-rem,]
standalone.TC.expr <- TC.expr[names(standalone.TCs),]

standalone.gene.TCs <- gene.TCs[intersect(names(standalone.TCs),names(gene.TCs)),]
standalone.gene.TC.expr <- gene.TC.expr[names(standalone.gene.TCs),]

## Chromatin state analysis
load("../data/ChIP_all_peaks.RData")
load("../data/ChIP_bivalent_peaks.RData")

## Assign chromatin states to TCs according to proximity with ChIP-seq peaks within 2kb
TC.states <- lapply(names(chip.k4.clone.peaks.grl), function(clone) {
    print(clone)
    state <- matrix(0,nrow=length(standalone.TCs),ncol=4)
    o <- findOverlaps(standalone.TCs,chip.k4.clone.peaks.grl[[clone]],maxgap=2000,type="any")
    state[unique(queryHits(o)),1] <- 1
    o <- findOverlaps(standalone.TCs,chip.k27.clone.peaks.grl[[clone]],maxgap=2000,type="any")
    state[unique(queryHits(o)),2] <- 1
    o <- findOverlaps(standalone.TCs,chip.k9.clone.peaks.grl[[clone]],maxgap=2000,type="any")
    state[unique(queryHits(o)),3] <- 1
    o <- findOverlaps(standalone.TCs,chip.bivalent.grl[[clone]],maxgap=2000,type="any")
    state[unique(queryHits(o)),1] <- 0
    state[unique(queryHits(o)),2] <- 0
    state[unique(queryHits(o)),4] <- 1

    state <- apply(state,1,function(x) paste(x,collapse=","))
    state[which(state=="0,0,0,0")] <- "none"
    state[which(state=="0,0,1,0")] <- "H3K9me3"
    state[which(state=="0,1,0,0")] <- "H3K27me3"
    state[which(state=="0,1,1,0")] <- "H3K27me3,H3K9me3"
    state[which(state=="1,0,0,0")] <- "H3K4me3"
    state[which(state=="1,0,1,0")] <- "H3K4me3,H3K9me3"
    state[which(state=="1,1,0,0")] <- "H3K4me3,H3K27me3"
    state[which(state=="1,1,1,0")] <- "H3K4me3,H3K27me3,H3K9me3"
    state[which(state=="0,0,0,1")] <- "H3K4me3+H3K27me3"
    state[which(state=="0,0,1,1")] <- "H3K4me3+H3K27me3,H3K9me3"
    names(state) <- names(standalone.TCs)
    state
})
names(TC.states) <- names(chip.k4.clone.peaks.grl)

## Focus TC sets
c410vs439.TCs <- rownames(clone.expressed)[which(clone.expressed[,"410"]>1 | clone.expressed[,"439"]>1)]
c410vs506vs438vs439.TCs <- rownames(clone.expressed)[which(apply(clone.expressed,1,function(x) any(x>1)))]

## Prepare frequency data frames

freq <- table(as.data.frame(cbind(TC.states[["410"]],TC.states[["439"]])[intersect(c410vs439.TCs,names(standalone.gene.TCs)),]))
df.c410vs439.gene <- data.frame(WT=rownames(freq),MUT=rep(colnames(freq),each=nrow(freq)),Freq=as.numeric(freq))

freq <- table(as.data.frame(cbind(TC.states[["410"]],TC.states[["439"]])[intersect(c410vs506vs438vs439.TCs,names(standalone.gene.TCs)),]))
df.c410vs439.gene.full <- data.frame(WT=rownames(freq),MUT=rep(colnames(freq),each=nrow(freq)),Freq=as.numeric(freq))

## Plot transition frequency (alluvial) plots

## Figure 1h
ggplot(df.c410vs439.gene, aes(y = Freq, axis1 = WT, axis2 = MUT)) +
    geom_alluvium(aes(fill = WT)) +
    geom_stratum(fill = "black", color = "grey") +
    geom_label(stat = "stratum", infer.label = TRUE) +
    scale_x_discrete(limits = c("WT (410)", "MCM2-2A (439)"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Paired") +
    ggtitle("Chromatin state at TSSs")
ggsave("plots/Chromatin_state_alluvial_410_vs_439_gene.pdf",width=8,height=16)

df.c410vs439.gene$fraction <- (df.c410vs439.gene$Freq+1) / sum(df.c410vs439.gene$Freq)

## Extended figure 3i
ggplot(df.c410vs439.gene,aes(x=MUT,y=WT,fill=fraction)) + geom_tile() + scale_fill_distiller(type = "seq", palette = "Blues", direction=1, trans = 'log10') +
    geom_text(aes(label=Freq),size=2) + xlab("439") + ylab("410") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/Chromatin_state_tile_410_vs_439_gene.pdf",width=6,height=5)

## Investigate chromatin state transitions with respect to RNA-seq differential expression status

load("../data/RNAseq_DE_genes.RData")

g.id <- as.character(sapply(gene.TCs$geneID,function(n) strsplit(n,"\\.")[[1]][1]))

down <- RNAseq.DE.genes.list[["full"]][["Down"]]
up <- RNAseq.DE.genes.list[["full"]][["Up"]]
down.TCs <- intersect(intersect(c410vs506vs438vs439.TCs,names(standalone.gene.TCs)),names(gene.TCs[which(g.id %in% down),]))
up.TCs <- intersect(intersect(c410vs506vs438vs439.TCs,names(standalone.gene.TCs)),names(gene.TCs[which(g.id %in% up),]))
all.TCs <- intersect(c410vs506vs438vs439.TCs,names(standalone.gene.TCs))
freq.down <- table(as.data.frame(cbind(TC.states[["410"]],TC.states[["439"]])[down.TCs,]))
freq.up <- table(as.data.frame(cbind(TC.states[["410"]],TC.states[["439"]])[up.TCs,]))

## Chromatin state transitions full model (WT vs MUT) RNA-seq DE

df.c410vs439.gene.DE.full <- rbind(data.frame(WT=rownames(freq.down),MUT=rep(colnames(freq.down),each=nrow(freq.down)),DE="Down",Freq=as.numeric(freq.down)),
                                   data.frame(WT=rownames(freq.up),MUT=rep(colnames(freq.up),each=nrow(freq.up)),DE="Up",Freq=as.numeric(freq.up)))

## Chromatin state enrichments

df.c410vs439.gene.DE.full$transition <- paste(as.character(df.c410vs439.gene.DE.full$WT),as.character(df.c410vs439.gene.DE.full$MUT),sep=" -> ")
df.c410vs439.gene.full$transition <- paste(as.character(df.c410vs439.gene.full$WT),as.character(df.c410vs439.gene.full$MUT),sep=" -> ")

u <- subset(df.c410vs439.gene.DE.full,DE=="Up")
d <- subset(df.c410vs439.gene.DE.full,DE=="Down")
a <- df.c410vs439.gene.full
wt.410vs439.full.state.enrich <- t(sapply(unique(df.c410vs439.gene.DE.full$WT), function(state) {
    print(state)
    um <- matrix(c(sum(subset(u,WT==state)$Freq),
                   sum(subset(u,WT!=state)$Freq),
                   sum(subset(a,WT==state)$Freq)-sum(subset(u,WT==state)$Freq),
                   sum(subset(a,WT!=state)$Freq)-sum(subset(u,WT!=state)$Freq)),
                 byrow=TRUE,ncol=2)
    ut <- fisher.test(um)
    dm <- matrix(c(sum(subset(d,WT==state)$Freq),
                   sum(subset(d,WT!=state)$Freq),
                   sum(subset(a,WT==state)$Freq)-sum(subset(d,WT==state)$Freq),
                   sum(subset(a,WT!=state)$Freq)-sum(subset(d,WT!=state)$Freq)),
                 byrow=TRUE,ncol=2)
    dt <- fisher.test(dm)
    c(ut$estimate,ut$p.value,dt$estimate,dt$p.value,sum(subset(u,WT==state)$Freq),sum(subset(d,WT==state)$Freq),sum(subset(a,WT==state)$Freq))
}))
colnames(wt.410vs439.full.state.enrich) <- c("Up, odds ratio", "Up, p-value", "Down, odds ratio", "Down, p-value","up","down","n")
rownames(wt.410vs439.full.state.enrich) <- unique(df.c410vs439.gene.DE.full$WT)

mut.410vs439.full.state.enrich <- t(sapply(unique(df.c410vs439.gene.DE.full$MUT), function(state) {
    um <- matrix(c(sum(subset(u,MUT==state)$Freq),
                   sum(subset(u,MUT!=state)$Freq),
                   sum(subset(a,MUT==state)$Freq)-sum(subset(u,MUT==state)$Freq),
                   sum(subset(a,MUT!=state)$Freq)-sum(subset(u,MUT!=state)$Freq)),
                 byrow=TRUE,ncol=2)
    ut <- fisher.test(um)
    dm <- matrix(c(sum(subset(d,MUT==state)$Freq),
                   sum(subset(d,MUT!=state)$Freq),
                   sum(subset(a,MUT==state)$Freq)-sum(subset(d,MUT==state)$Freq),
                   sum(subset(a,MUT!=state)$Freq)-sum(subset(d,MUT!=state)$Freq)),
                 byrow=TRUE,ncol=2)
    dt <- fisher.test(dm)
    c(ut$estimate,ut$p.value,dt$estimate,dt$p.value,sum(subset(u,WT==state)$Freq),sum(subset(d,WT==state)$Freq),sum(subset(a,WT==state)$Freq))
}))
colnames(mut.410vs439.full.state.enrich) <- c("Up, odds ratio", "Up, p-value", "Down, odds ratio", "Down, p-value","up","down","n")
rownames(mut.410vs439.full.state.enrich) <- unique(df.c410vs439.gene.DE.full$MUT)

c410vs439.full.transition.enrich <- t(sapply(unique(df.c410vs439.gene.DE.full$transition), function(t) {
    um <- matrix(c(sum(subset(u,transition==t)$Freq),
                   sum(subset(u,transition!=t)$Freq),
                   sum(subset(a,transition==t)$Freq)-sum(subset(u,transition==t)$Freq),
                   sum(subset(a,transition!=t)$Freq)-sum(subset(u,transition!=t)$Freq)),
                 byrow=TRUE,ncol=2)
    ut <- fisher.test(um)
    dm <- matrix(c(sum(subset(d,transition==t)$Freq),
                   sum(subset(d,transition!=t)$Freq),
                   sum(subset(a,transition==t)$Freq)-sum(subset(d,transition==t)$Freq),
                   sum(subset(a,transition!=t)$Freq)-sum(subset(d,transition!=t)$Freq)),
                 byrow=TRUE,ncol=2)
    dt <- fisher.test(dm)
    c(ut$estimate,ut$p.value,dt$estimate,dt$p.value,sum(subset(u,transition==t)$Freq),sum(subset(d,transition==t)$Freq),sum(subset(a,transition==t)$Freq))
}))
colnames(c410vs439.full.transition.enrich) <- c("Up, odds ratio", "Up, p-value", "Down, odds ratio", "Down, p-value","up","down","n")
rownames(c410vs439.full.transition.enrich) <- unique(df.c410vs439.gene.DE.full$transition)

df <- rbind(data.frame(p.value=c(mut.410vs439.full.state.enrich[,"Up, p-value"],mut.410vs439.full.state.enrich[,"Down, p-value"]),
                       odds.ratio=c(mut.410vs439.full.state.enrich[,"Up, odds ratio"],mut.410vs439.full.state.enrich[,"Down, odds ratio"]),
                       count=c(mut.410vs439.full.state.enrich[,"up"],mut.410vs439.full.state.enrich[,"down"]),
                       state=rownames(mut.410vs439.full.state.enrich),
                       DE=rep(c("Up","Down"),each=nrow(mut.410vs439.full.state.enrich)),
                       clone="439"),
            data.frame(p.value=c(wt.410vs439.full.state.enrich[,"Up, p-value"],wt.410vs439.full.state.enrich[,"Down, p-value"]),
                       odds.ratio=c(wt.410vs439.full.state.enrich[,"Up, odds ratio"],wt.410vs439.full.state.enrich[,"Down, odds ratio"]),
                       count=c(wt.410vs439.full.state.enrich[,"up"],wt.410vs439.full.state.enrich[,"down"]),
                       state=rownames(wt.410vs439.full.state.enrich),
                       DE=rep(c("Up","Down"),each=nrow(wt.410vs439.full.state.enrich)),
                       clone="410"))

df[df$p.value>0.05,"odds.ratio"] <- NA

## Figure 2c
ggplot(subset(df,clone="410"),aes(x=DE,y=state,fill=log2(odds.ratio))) + geom_tile() + scale_fill_distiller(type = "div", palette = "RdBu") +
    geom_text(aes(label=count),size=2) + xlab("") + ylab("")
ggsave("plots/410_vs_439_WT_chromatin_state_enrich_full_DE_model_tile.pdf",width=4,height=4)

df <- data.frame(p.value=c(c410vs439.full.transition.enrich[,"Up, p-value"],c410vs439.full.transition.enrich[,"Down, p-value"]),
                 odds.ratio=c(c410vs439.full.transition.enrich[,"Up, odds ratio"],c410vs439.full.transition.enrich[,"Down, odds ratio"]),
                 count=c(c410vs439.full.transition.enrich[,"up"],c410vs439.full.transition.enrich[,"down"]),
                 transition=rownames(c410vs439.full.transition.enrich),
                 WT=as.character(sapply(rownames(c410vs439.full.transition.enrich),function(n) strsplit(n," -> ",)[[1]][1])),
                 MUT=as.character(sapply(rownames(c410vs439.full.transition.enrich),function(n) strsplit(n," -> ",)[[1]][2])),
                 n=c410vs439.full.transition.enrich[,"n"],
                 DE=rep(c("Up","Down"),each=nrow(c410vs439.full.transition.enrich)))

tmp <- subset(df,p.value<=0.05 & count>=10)
keep <- unique(c(as.character(tmp$WT),as.character(tmp$MUT)))
df <- subset(df,WT %in% keep)
df <- subset(df,MUT %in% keep)
df[df$p.value>0.05 | df$count < 10,"odds.ratio"]  <- NA

## Extended figure 5b
ggplot(df,aes(x=MUT,y=WT,fill=log2(odds.ratio))) + geom_tile() + scale_fill_distiller(type = "div", palette = "RdBu") +
    geom_text(aes(label=count),size=2) + xlab("") + ylab("") + facet_wrap(~DE,ncol=2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/410_vs_439_transition_chromatin_state_enrich_full_DE_model_tile.pdf",width=9,height=5)


## ChIP-seq signal at bivalent peaks vs DE status
load("../data//ChIP_counts_wt_peaks.RData")

o <- findOverlaps(standalone.TCs,chip.bivalent.grl[["410"]],maxgap=2000,type="any")
df.TCs.biv.peaks.410.439.gene.full <- data.frame(TC=names(standalone.TCs)[queryHits(o)],
                                                 peak=names(chip.bivalent.grl[["410"]])[subjectHits(o)],
                                                 WT.state=TC.states[["410"]][names(standalone.TCs)[queryHits(o)]],
                                                 MUT.state=TC.states[["439"]][names(standalone.TCs)[queryHits(o)]],
                                                 stringsAsFactors=FALSE)

down <- RNAseq.DE.genes.list[["full"]][["Down"]]
up <- RNAseq.DE.genes.list[["full"]][["Up"]]
down.TCs <- intersect(intersect(c410vs506vs438vs439.TCs,names(standalone.gene.TCs)),names(gene.TCs[which(g.id %in% down),]))
up.TCs <- intersect(intersect(c410vs506vs438vs439.TCs,names(standalone.gene.TCs)),names(gene.TCs[which(g.id %in% up),]))
all.TCs <- intersect(c410vs506vs438vs439.TCs,names(standalone.gene.TCs))

df.TCs.biv.peaks.410.439.gene.full <- subset(df.TCs.biv.peaks.410.439.gene.full,TC %in% all.TCs)

df.TCs.biv.peaks.410.439.gene.full$DE <- "constant"
df.TCs.biv.peaks.410.439.gene.full[df.TCs.biv.peaks.410.439.gene.full$TC %in% down.TCs,"DE"] <- "Down"
df.TCs.biv.peaks.410.439.gene.full[df.TCs.biv.peaks.410.439.gene.full$TC %in% up.TCs,"DE"] <- "Up"

focus.data <- assay(chip.peak.se.list$K4.410,"counts")[unique(df.TCs.biv.peaks.410.439.gene.full$peak),]
dimnames(focus.data) <- list(unique(df.TCs.biv.peaks.410.439.gene.full$peak),colnames(chip.peak.se.list$K4.410))
focus.data.RRPM <- t(apply(focus.data,1,function(x) x*colData(chip.peak.se.list$K4.410)$RRPM.scaling))

data.type <- as.character(sapply(colnames(focus.data),function(n) strsplit(n,"_ChIP")[[1]][1]))
focus.data.mean.RRPM <- sapply(unique(data.type),function(n) rowMeans(focus.data.RRPM[,data.type==n]))
dimnames(focus.data.mean.RRPM) <- list(rownames(focus.data.RRPM),unique(data.type))

df.TCs.biv.peaks.410.439.gene.full.ChIP.ratio  <- rbind(df.TCs.biv.peaks.410.439.gene.full,df.TCs.biv.peaks.410.439.gene.full,df.TCs.biv.peaks.410.439.gene.full)
df.TCs.biv.peaks.410.439.gene.full.ChIP.ratio$FC <- c(log2(focus.data.mean.RRPM[as.character(df.TCs.biv.peaks.410.439.gene.full$peak),"439_H3K4me3"]+1) -
                                                      log2(focus.data.mean.RRPM[as.character(df.TCs.biv.peaks.410.439.gene.full$peak),"410_H3K4me3"]+1),
                                                      log2(focus.data.mean.RRPM[as.character(df.TCs.biv.peaks.410.439.gene.full$peak),"439_H3K27me3"]+1) -
                                                      log2(focus.data.mean.RRPM[as.character(df.TCs.biv.peaks.410.439.gene.full$peak),"410_H3K27me3"]+1),
                                                      log2(focus.data.mean.RRPM[as.character(df.TCs.biv.peaks.410.439.gene.full$peak),"439_H3K9me3"]+1) -
                                                      log2(focus.data.mean.RRPM[as.character(df.TCs.biv.peaks.410.439.gene.full$peak),"410_H3K9me3"]+1))

df.TCs.biv.peaks.410.439.gene.full.ChIP.ratio$ChIP <- rep(c("H3K4me3","H3K27me3","H3K9me3"),each=nrow(df.TCs.biv.peaks.410.439.gene.full))

## Figure 2d
ggplot(subset(df.TCs.biv.peaks.410.439.gene.full.ChIP.ratio,WT.state=="H3K4me3+H3K27me3" & MUT.state=="H3K4me3+H3K27me3" & ChIP %in% c("H3K4me3","H3K27me3") & DE %in% c("Up","Down")),
       aes(ChIP,FC,fill=DE)) + geom_boxplot() +
    xlab("") + ylab("log2 fold change 439 vs 410")
ggsave("plots/410_bivalent_ChIP_ratio_439_gene_DE_full_subset.pdf",width=4,height=4)
