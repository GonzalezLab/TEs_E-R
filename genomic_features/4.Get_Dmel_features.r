#!/usr/bin/env Rscript

library(GenomicFeatures)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with the dmel-all-r6.35.gtf file downloaded from Flybase (i.e., "/home/")
#args[2]: path to the results folder

parent.f <- args[1]
parent.o <- args[2]

# Read gft file
txdb <- makeTxDbFromGFF(paste(parent.f, "dmel-all-r6.35.gtf", sep=""))

# Subset for chromosomes
seqlevels(txdb) <- c("2L","2R","3L","3R","4","X","Y")
 
# GRanges for genes
genes <- genes(txdb)

# GRanges for each exon of each gene
ex <- exonsBy(txdb,by="gene") 

# GRanges for each intron of each gene
int <- intronsByTranscript(txdb, use.names=TRUE) 

# GRanges for each 1000 bp upstream the transcription start site
pro <- promoters(txdb, upstream=1000, downstream=0)

# GRanges for 5'UTR and 3'UTR of each gene
utr3 <- threeUTRsByTranscript(txdb,use.names=TRUE)
utr5 <- fiveUTRsByTranscript(txdb,use.names=TRUE)

# Remove redundancy in ranges
ex.1 <- as.data.frame(ex)
ex.2 <- GRanges(ex.1$seqnames, IRanges(ex.1$start, ex.1$end))
ex.3 <- reduce(ex.2)
ex.f <- as.data.frame(ex.3)

int.1 <- as.data.frame(int)
int.2 <- GRanges(int.1$seqnames, IRanges(int.1$start, int.1$end))
int.3 <- reduce(int.2)
int.f <- as.data.frame(int.3)

pro.1 <- as.data.frame(pro)
pro.2 <- GRanges(pro.1$seqnames, IRanges(pro.1$start, pro.1$end))
pro.3 <- reduce(pro.2)
pro.f <- as.data.frame(pro.3)

utr5.1 <- as.data.frame(utr5)
utr5.2 <- GRanges(utr5.1$seqnames, IRanges(utr5.1$start, utr5.1$end))
utr5.3 <- reduce(utr5.2)
utr5.f <- as.data.frame(utr5.3)

utr3.1 <- as.data.frame(utr3)
utr3.2 <- GRanges(utr3.1$seqnames, IRanges(utr3.1$start, utr3.1$end))
utr3.3 <- reduce(utr3.2)
utr3.f <- as.data.frame(utr3.3)

# Get 1000 bp from 5' and 3' of genes
genes.c <- data.frame(genes)

u5 <- u3 <- u5.all <- u3.all <- NULL
for (i in 1:length(genes.c[,1])) {
	if (genes.c[i, 5] == "+") {
		u5 <- cbind(as.character(genes.c[i,1]), genes.c[i,2] - 1000, genes.c[i,2], genes.c[i,4], "5UPSTR")
		u3 <- cbind(as.character(genes.c[i,1]), genes.c[i,3], genes.c[i,3] + 1000, genes.c[i,4], "3UPSTR")
		u5.all <- rbind(u5.all, u5)
		u3.all <- rbind(u3.all, u3)
	}
	if (genes.c[i, 5] == "-") {
		u5 <- cbind(as.character(genes.c[i,1]), genes.c[i,3], genes.c[i,3] + 1000, genes.c[i,4], "5UPSTR")
		u3 <- cbind(as.character(genes.c[i,1]), genes.c[i,2] - 1000, genes.c[i,2], genes.c[i,4], "3UPSTR")
		u5.all <- rbind(u5.all, u5)
		u3.all <- rbind(u3.all, u3)
	}	
}
u5.all <- data.frame(u5.all,stringsAsFactors=F)
u3.all <- data.frame(u3.all,stringsAsFactors=F)	
colnames(u5.all) <- c("chr","start","end","width","features")
colnames(u3.all) <- c("chr","start","end","width","features")

# Merge all (with 5upstr and 3upstr)
ex.f.1 <- cbind(ex.f[,1:4], rep("exon", dim(ex.f)[1]))
int.f.1 <- cbind(int.f[,1:4], rep("intron", dim(int.f)[1]))
pro.f.1 <- cbind(pro.f[,1:4], rep("promoter", dim(pro.f)[1]))
utr5.f.1 <- cbind(utr5.f[,1:4], rep("5'UTR", dim(utr5.f)[1]))
utr3.f.1 <- cbind(utr3.f[,1:4], rep("3'UTR", dim(utr3.f)[1]))

colnames(ex.f.1) <- c("chr","start","end","width","features")
colnames(int.f.1) <- c("chr","start","end","width","features")
colnames(pro.f.1) <- c("chr","start","end","width","features")
colnames(utr5.f.1) <- c("chr","start","end","width","features")
colnames(utr3.f.1) <- c("chr","start","end","width","features")

all <- rbind(ex.f.1,int.f.1,pro.f.1,utr5.f.1,utr3.f.1,u5.all,u3.all)   
all.f <- all[order(as.character(all$chr), all$start, all$end), ]

# Change the name of chromosomes
all.f$chr <- gsub("2L", "chr2L", all.f$chr)
all.f$chr <- gsub("2R", "chr2R", all.f$chr)
all.f$chr <- gsub("3L", "chr3L", all.f$chr)
all.f$chr <- gsub("3R", "chr3R", all.f$chr)
all.f$chr <- gsub("4", "chr4", all.f$chr)
all.f$chr <- gsub("X", "chrX", all.f$chr)
all.f$chr <- gsub("Y", "chrY", all.f$chr)

write.table(all.f, file=paste(parent.o, "dmel-all-r6.35_mod2.gtf", sep=""), quote=F, sep="\t", col.names=T, row.names=F)

