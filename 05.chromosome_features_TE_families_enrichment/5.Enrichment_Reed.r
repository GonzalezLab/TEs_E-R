#!/usr/bin/env Rscript

library(DescTools) 
library(GenomicRanges)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with the diferent subfolders (pools) containing the TEMP output files (i.e., "/home/")
#args[2]: path to the results folder
#args[3]: file with the TEs dictionary FlyBase-Repbase (with full path)
#args[4]: gtf file ("dmel-all-r6.35_mod2.gtf", with full path)

parent.f <- args[1]
parent.o <- args[2]
dic <- args[3]
gtf <- args[4]

# Read dictionary
d <- read.table(file=dic, header=T, stringsAsFactors=F)

## Read files with predicted TEs (the columns are: Chr, TE_ID, Start_BaseLine, Stop_BaseLine, Start_Selection_pool, Stop_Selection_pool, Frequency_BaseLine, Frequency_Selection_pool, Pval, Pools_compared)
# TEs predicted with TEMP and TIDAL
pool1_tt <- read.table(file=paste(parent.f, "reed_intervals_concatenated_1p1.f0.10.txt", sep=""), 
						header=T, stringsAsFactors=F)

# TEs predicted with only TEMP
pool1_t <- read.table(file=paste(parent.f, "only_temp/", "reed_intervals_concatenated_1p1.f0.10_temp.txt", 
						sep=""), header=T, stringsAsFactors=F)

## Remove significant TEs
# Significant TEs predicted with TEMP and TIDAL (the columns are: chr, TE_ID, start, stop, number of the pools where the TE is significant)
sig_tt <- sig_tt <- read.table(file=paste(parent.f, "significant_reed.txt", sep=""), 
								header=F, stringsAsFactors=F)

ind <- ind2 <- NULL
for (i in 1:dim(sig_tt)[1]) {
	ind <- which(sig_tt[i,1] == pool1_tt[,1] & sig_tt[i,2] == pool1_tt[,2])
	for (j in 1:length(ind)) {
	 	if (c(pool1_tt[ind[j],3], pool1_tt[ind[j],4]) %overlaps% c(sig_tt[i,3], sig_tt[i,4]) == TRUE) {
	 		ind2 <- c(ind2, ind[j])
	 	}	
	 }		
}

if (length(ind2) == dim(sig_tt)[1]) {
	pool1_tt <- pool1_tt[-ind2, ]	
} else { 
	ind3 <- which(pool1_tt[ind2, 9] < 0.05)
	pool1_tt <- pool1_tt[-ind2[ind3], ]
}
	
# Significant TEs predicted with only TEMP
sig_t <- read.table(file=paste(parent.f, Study,"only_temp/", "significant_reed_temp.txt", sep=""), 
					header=F, stringsAsFactors=F)
ind <- ind2 <- NULL
for (i in 1:dim(sig_t)[1]) {
	ind <- which(sig_t[i,1] == pool1_t[,1] & sig_t[i,2] == pool1_t[,2])
	for (j in 1:length(ind)) {
	 	if (c(pool1_t[ind[j],3], pool1_t[ind[j],4]) %overlaps% c(sig_t[i,3], sig_t[i,4]) == TRUE) {
	 		ind2 <- c(ind2, ind[j])
	 	}	
	 }		
}

if (length(ind2) == dim(sig_t)[1]) {
	pool1_t <- pool1_t[-ind2, ]	
} else { 
	ind3 <- which(pool1_t[ind2, 9] < 0.05)
	pool1_t <- pool1_t[-ind2[ind3], ]
}

## Replace the name of the TE using the dictionary
# One pool
for (i in 1:dim(pool1_t)[1]) {
	if (pool1_t[i,2] %in% d[,2]) {
		pool1_t[i,2] <- as.character(d[d[,2] == pool1_t[i,2], 1])
	}	
}	
	
# Find unique rows
pool1_tt <- unique(pool1_tt[, c("Chr","TE_ID","Start_BL","Stop_BL","Start_S","Stop_S","Freq_BL",
								"Freq_S","Pval","Comparison")])
pool1_t <- unique(pool1_t[, c("Chr","TE_ID","Start_BL","Stop_BL","Start_S","Stop_S","Freq_BL",
								"Freq_S","Pval","Comparison")])

# Remove redundancy between those detected by two methods and those detected only with temp
df1 <- sm1 <- NULL
for (i in 1:length(pool1_tt[,1])) {
	c <- 0
	for (j in 1:length(pool1_t[,1])) {
		if ( pool1_tt[i,1] == pool1_t[j,1] & pool1_tt[i,2] == pool1_t[j,2] ) {
			if (c(pool1_tt[i,5], pool1_tt[i,6]) %overlaps% c(pool1_t[j,5], pool1_t[j,6]) == TRUE) {	
				c <- c + 1
			}
		}
	}
	if (c == 0) {
		df1 <- rbind(df1, pool1_tt[i, ])
	} else { 
		sm1 <- rbind(sm1, pool1_tt[i, ]) 
	}
}

df2 <- sm2 <- NULL
for (i in 1:length(pool1_t[,1])) {
	c <- 0
	for (j in 1:length(pool1_tt[,1])) {
		if ( pool1_t[i,1] == pool1_tt[j,1] & pool1_t[i,2] == pool1_tt[j,2] ) {
			if (c(pool1_t[i,5], pool1_t[i,6]) %overlaps% c(pool1_tt[j,5], pool1_tt[j,6]) == TRUE) {	
				c <- c + 1
			}
		}
	}
	if (c == 0) {
		df2 <- rbind(df2, pool1_t[i, ])
	} else { 
		sm2 <- rbind(sm2, pool1_t[i, ]) 
	}
}			

# Merge the rows that are different between TEMP-TIDAL and TEMP-only files and the rows that have overlapping TE coordinates
all <- NULL
if (dim(sm1)[1] == dim(sm2)[1] | dim(sm1)[1] > dim(sm2)[1]) {
	all <- rbind(df1, sm1, df2)
} else {
 	all <- rbind(df1, sm2, df2)
}

## Read files with significant TEs
# Significant TEs predicted with TEMP and TIDAL
sig_tt <- read.table(file=paste(parent.f, "significant_reed.txt", sep=""), header=F, stringsAsFactors=F)
sig1_tt <- sig_tt[sig_tt[ ,5 ] == 1, ]

# Significant TEs predicted with only TEMP
sig_t <- read.table(file=paste(parent.f, "only_temp/", "significant_reed_temp.txt", sep=""), 
					header=F, stringsAsFactors=F)

# Replace the name of the TE using the dictionary
for (i in 1:dim(sig_t)[1]) {
	if (sig_t[i,2] %in% d[,2]) {
		sig_t[i,2] <- as.character(d[d[,2] == sig_t[i,2], 1])
	}	
}
  
sig1_t <- sig_t[sig_t[ ,5] == 1, ]

# Find unique rows
sig1_tt <- unique(sig1_tt[, c("V1","V2","V3","V4","V5")])
sig1_t <- unique(sig1_t[, c("V1","V2","V3","V4","V5")])

# Check if significant TEs in TEMP-TIDAL one pool are detected by TEMP-only in two-three pools
# File with TIDAL-TEMP one pool as query and TEMP-only more than one pool as subject
sdf1 <- ssm1 <- NULL
for (i in 1:length(sig1_tt[,1])) {
	c <- 0
	for (j in 1:length(sig_t[,1])) {
		if ( sig1_tt[i,1] == sig_t[j,1] & sig1_tt[i,2] == sig_t[j,2] ) {
			if (c(sig1_tt[i,3], sig1_tt[i,4]) %overlaps% c(sig_t[j,3], sig_t[j,4]) == TRUE) {	
				c <- c + 1
			}
		}
	}
	if (c == 0) {
		sdf1 <- rbind(sdf1, sig1_tt[i, ])
	} else { 
		ssm1 <- rbind(ssm1, sig1_tt[i, ]) 
	}
}			

# TEMP-only more than one pool as query and TIDAL-TEMP one pool as subject
sdf2 <- ssm2 <- NULL
for (i in 1:length(sig_t[,1])) {
	c <- 0
	for (j in 1:length(sig1_tt[,1])) {
		if ( sig_t[i,1] == sig1_tt[j,1] & sig_t[i,2] == sig1_tt[j,2] ) {
			if (c(sig_t[i,3], sig_t[i,4]) %overlaps% c(sig1_tt[j,3], sig1_tt[j,4]) == TRUE) {	
				c <- c + 1
			}
		}
	}
	if (c == 0) {
		sdf2 <- rbind(sdf2, sig_t[i, ])
	} else { 
		ssm2 <- rbind(ssm2, sig_t[i, ]) 
	}
}	

# Merge the rows that are different between TEMP-TIDAL and TEMP-only files and the rows that have overlapping TE coordinates
all.s <- NULL
if (dim(unique(ssm1))[1] == dim(unique(ssm1))[1] | dim(unique(ssm1))[1] > dim(unique(ssm2))[1]) {
	all.s <- rbind(sdf1, ssm1, sdf2)
} else {
 	all.s <- rbind(sdf1, ssm2, sdf2)
}

# Merge all the SE and WE significant and non-significant TEs
we_all <- all.s
ns_we_all <- all

# Read gtf file (the columns are: chr, start, stop, width, features)
gff <- read.csv(file=gtf, header=T, sep="\t", stringsAsFactors=F)
we_all.f <- ns_we_all.f <- NULL
	
# Get the chromosome features where the significant WE TEs where located
we_all.f <- cbind(we_all, rep(NA, length(we_all[,1])))
for (i in 1:length(we_all[,1])) {
	sub.gff <- gff[gff[,1] == we_all[i, 1], ]
	o <- NULL
	for (j in 1:length(sub.gff[,1])) {
		if (c(as.numeric(we_all[i,3]), as.numeric(we_all[i,4])) %overlaps% c(sub.gff[j,2], sub.gff[j,3]) == TRUE) {	
			o <- rbind(o, sub.gff[j, ])
		}
	}
	if (length(unique(o[,5])) > 1 & (is.element("5UPSTR", unique(o[,5])) | is.element("3UPSTR", unique(o[,5])))) {
		ft <- unique(o[,5])[!unique(o[,5]) %in% "3UPSTR" & !unique(o[,5]) %in% "5UPSTR"]
	} else {
		ft <- unique(o[,5])	
	}	
	we_all.f[i, 6] <- paste(ft, collapse=';')
}

colnames(we_all.f) <- c("Chr", "TE_ID", "Start", "Stop", "No pools", "features")
write.table(we_all.f, file=paste(parent.o, "significant_we_tes_reed_features_gtf2.txt", sep=""), 
			sep="\t", quote=F, col.names=T, row.names=F)

# Read the file
we_all.f <- read.table(file=paste(parent.o, "significant_we_tes_reed_features_gtf2.txt", sep=""), 
						header=T, quote="", sep="\t", stringsAsFactors=F)
	
# Get chromosome features where the non-significant WE TEs where located
ns_we_all.f <- cbind(ns_we_all, rep(NA, length(ns_we_all[,1])))

for (i in 1:length(ns_we_all[,1])) {
	sub.gff <- gff[gff[,1] == ns_we_all[i, 1], ]
	o <- NULL
	for (j in 1:length(sub.gff[,1])) {
		if (c(as.numeric(ns_we_all[i,5]), as.numeric(ns_we_all[i,6])) %overlaps% c(sub.gff[j,2], sub.gff[j,3]) == TRUE) {	
			o <- rbind(o,sub.gff[j, ])
		}
	}
	if (length(unique(o[,5])) > 1 & (is.element("5UPSTR", unique(o[,5])) | is.element("3UPSTR", unique(o[,5])))) {
		ft <- unique(o[,5])[!unique(o[,5]) %in% "3UPSTR" & !unique(o[,5]) %in% "5UPSTR"]
	} else {
		ft <- unique(o[,5])	
	}	
	ns_we_all.f[i, 5] <- paste(ft, collapse=';')
}	
colnames(ns_we_all.f) <- c("Chr", "TE_ID", "Start", "Stop", "features") 
write.table(ns_we_all.f, file=paste(parent.o, "non-significant_we_tes_reed_features_gtf2.txt", sep=""), 
			sep="\t", quote=F, col.names=T, row.names=F)	

#Read the file
ns_we_all.f <- read.table(file=paste(parent.o, "non-significant_we_tes_reed_features_gtf2.txt", sep=""), 
							quote="", header=T, sep="\t", stringsAsFactors=F)

#### Perform enrichment analysis (chromosome and features) for WE TEs
# Count the number of significant and non-significant TEs in each feature
we_all.f.f <- length(grep("5'UTR",we_all.f[,6]))
ns_we_all.f.f <- length(grep("5'UTR",ns_we_all.f[,5]))
we_all.f.nf <- dim(we_all.f)[1] - length(grep("5'UTR",we_all.f[,6]))
ns_we_all.f.nf <- dim(ns_we_all.f)[1] - length(grep("5'UTR",ns_we_all.f[,5]))

we_all.f.e <- length(grep("exon",we_all.f[,6]))
ns_we_all.f.e <- length(grep("exon",ns_we_all.f[,5]))
we_all.f.ne <- dim(we_all.f)[1] - length(grep("exon",we_all.f[,6]))
ns_we_all.f.ne <- dim(ns_we_all.f)[1] - length(grep("exon",ns_we_all.f[,5]))

we_all.f.i <- length(grep("intron",we_all.f[,6]))
ns_we_all.f.i <- length(grep("intron",ns_we_all.f[,5]))
we_all.f.ni <- dim(we_all.f)[1] - length(grep("intron",we_all.f[,6]))
ns_we_all.f.ni <- dim(ns_we_all.f)[1] - length(grep("intron",ns_we_all.f[,5]))

we_all.f.t <- length(grep("3'UTR",we_all.f[,6]))
ns_we_all.f.t <- length(grep("3'UTR",ns_we_all.f[,5]))
we_all.f.nt <- dim(we_all.f)[1] - length(grep("3'UTR",we_all.f[,6]))
ns_we_all.f.nt <- dim(ns_we_all.f)[1] - length(grep("3'UTR",ns_we_all.f[,5]))

we_all.f.r <- length(grep("promoter",we_all.f[,6]))
ns_we_all.f.r <- length(grep("promoter",ns_we_all.f[,5]))
we_all.f.nr <- dim(we_all.f)[1] - length(grep("promoter",we_all.f[,6]))
ns_we_all.f.nr <- dim(ns_we_all.f)[1] - length(grep("promoter",ns_we_all.f[,5]))

we_all.f.ur <- length(grep("5UPSTR",we_all.f[,6]))
ns_we_all.f.ur <- length(grep("5UPSTR",ns_we_all.f[,5]))
we_all.f.nur <- dim(we_all.f)[1] - length(grep("5UPSTR",we_all.f[,6]))
ns_we_all.f.nur <- dim(ns_we_all.f)[1] - length(grep("5UPSTR",ns_we_all.f[,5]))

we_all.f.dr <- length(grep("3UPSTR",we_all.f[,6]))
ns_we_all.f.dr <- length(grep("3UPSTR",ns_we_all.f[,5]))
we_all.f.ndr <- dim(we_all.f)[1] - length(grep("3UPSTR",we_all.f[,6]))
ns_we_all.f.ndr <- dim(ns_we_all.f)[1] - length(grep("3UPSTR",ns_we_all.f[,5]))

# Count the number of significant and non-significant WE TEs in each chromosomal arm
we_all.f.2l <- length(grep("^chr2L",we_all.f[,1]))
ns_we_all.f.2l <- length(grep("^chr2L",ns_we_all.f[,1]))
we_all.f.n2l <- dim(we_all.f)[1] - length(grep("^chr2L",we_all.f[,1]))
ns_we_all.f.n2l <- dim(ns_we_all.f)[1] - length(grep("^chr2L",ns_we_all.f[,1]))

we_all.f.2r <- length(grep("^chr2R",we_all.f[,1]))
ns_we_all.f.2r <- length(grep("^chr2R",ns_we_all.f[,1]))
we_all.f.n2r <- dim(we_all.f)[1] - length(grep("^chr2R",we_all.f[,1]))
ns_we_all.f.n2r <- dim(ns_we_all.f)[1] - length(grep("^chr2R",ns_we_all.f[,1]))

we_all.f.3l <- length(grep("^chr3L",we_all.f[,1]))
ns_we_all.f.3l <- length(grep("^chr3L",ns_we_all.f[,1]))
we_all.f.n3l <- dim(we_all.f)[1] - length(grep("^chr3L",we_all.f[,1]))
ns_we_all.f.n3l <- dim(ns_we_all.f)[1] - length(grep("^chr3L",ns_we_all.f[,1]))

we_all.f.3r <- length(grep("^chr3R",we_all.f[,1]))
ns_we_all.f.3r <- length(grep("^chr3R",ns_we_all.f[,1]))
we_all.f.n3r <- dim(we_all.f)[1] - length(grep("^chr3R",we_all.f[,1]))
ns_we_all.f.n3r <- dim(ns_we_all.f)[1] - length(grep("^chr3R",ns_we_all.f[,1]))

we_all.f.4 <- length(grep("^chr4",we_all.f[,1]))
ns_we_all.f.4 <- length(grep("^chr4",ns_we_all.f[,1]))
we_all.f.n4 <- dim(we_all.f)[1] - length(grep("^chr4",we_all.f[,1]))
ns_we_all.f.n4 <- dim(ns_we_all.f)[1] - length(grep("^chr4",ns_we_all.f[,1]))

we_all.f.x <- length(grep("^chrX",we_all.f[,1]))
ns_we_all.f.x <- length(grep("^chrX",ns_we_all.f[,1]))
we_all.f.nx <- dim(we_all.f)[1] - length(grep("^chrX",we_all.f[,1]))
ns_we_all.f.nx <- dim(ns_we_all.f)[1] - length(grep("^chrX",ns_we_all.f[,1]))

we_all.f.y <- length(grep("^chrY",we_all.f[,1]))
ns_we_all.f.y <- length(grep("^chrY",we_all.f[,1]))
we_all.f.ny <- dim(we_all.f)[1] - length(grep("^chrY",we_all.f[,1]))
ns_we_all.f.ny <- dim(we_all.f)[1] - length(grep("^chrY",we_all.f[,1]))

f.fisher <- fisher.test(matrix(c(we_all.f.f, ns_we_all.f.f, we_all.f.nf, ns_we_all.f.nf), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
e.fisher <- fisher.test(matrix(c(we_all.f.e, ns_we_all.f.e, we_all.f.ne, ns_we_all.f.ne), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
i.fisher <- fisher.test(matrix(c(we_all.f.i, ns_we_all.f.i, we_all.f.ni, ns_we_all.f.ni), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
t.fisher <- fisher.test(matrix(c(we_all.f.t, ns_we_all.f.t, we_all.f.nt, ns_we_all.f.nt), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
r.fisher <- fisher.test(matrix(c(we_all.f.r, ns_we_all.f.r, we_all.f.nr, ns_we_all.f.nr), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
ur.fisher <- fisher.test(matrix(c(we_all.f.ur, ns_we_all.f.ur, we_all.f.nur, ns_we_all.f.nur), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
dr.fisher <- fisher.test(matrix(c(we_all.f.dr, ns_we_all.f.dr, we_all.f.ndr, ns_we_all.f.ndr), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value

fisher.2l <- fisher.test(matrix(c(we_all.f.2l, ns_we_all.f.2l, we_all.f.n2l, ns_we_all.f.n2l), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
fisher.2r <- fisher.test(matrix(c(we_all.f.2r, ns_we_all.f.2r, we_all.f.n2r, ns_we_all.f.n2r), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
fisher.3l <- fisher.test(matrix(c(we_all.f.3l, ns_we_all.f.3l, we_all.f.n3l, ns_we_all.f.n3l), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
fisher.3r <- fisher.test(matrix(c(we_all.f.3r, ns_we_all.f.3r, we_all.f.n3r, ns_we_all.f.n3r), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
fisher.4 <- fisher.test(matrix(c(we_all.f.4, ns_we_all.f.4, we_all.f.n4, ns_we_all.f.n4), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
fisher.x <- fisher.test(matrix(c(we_all.f.x, ns_we_all.f.x, we_all.f.nx, ns_we_all.f.nx), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value
fisher.y <- fisher.test(matrix(c(we_all.f.y, ns_we_all.f.y, we_all.f.ny, ns_we_all.f.ny), ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value

pvals <- c(f.fisher, e.fisher, i.fisher, t.fisher, r.fisher, ur.fisher, dr.fisher,
	fisher.2l,fisher.2r,fisher.3l,fisher.3r,fisher.4,fisher.x,fisher.y)

qval <- p.adjust(pvals, method = "fdr", n = length(pvals)) 
counts <- rbind(cbind(ns_we_all.f.nf, we_all.f.nf, ns_we_all.f.f, we_all.f.f),
				cbind(ns_we_all.f.ne, we_all.f.ne, ns_we_all.f.e, we_all.f.e),
				cbind(ns_we_all.f.ni, we_all.f.ni, ns_we_all.f.i, we_all.f.i),
				cbind(ns_we_all.f.nt, we_all.f.nt, ns_we_all.f.t, we_all.f.t),
				cbind(ns_we_all.f.nr, we_all.f.nr, ns_we_all.f.r, we_all.f.r),
 				cbind(ns_we_all.f.nur, we_all.f.nur, ns_we_all.f.ur, we_all.f.ur),
 				cbind(ns_we_all.f.ndr, we_all.f.ndr, ns_we_all.f.dr, we_all.f.dr),
 				cbind(ns_we_all.f.n2l, we_all.f.n2l, ns_we_all.f.2l, we_all.f.2l),
 				cbind(ns_we_all.f.n2r, we_all.f.n2r, ns_we_all.f.2r, we_all.f.2r),
 				cbind(ns_we_all.f.n3l, we_all.f.n3l, ns_we_all.f.3l, we_all.f.3l),
 				cbind(ns_we_all.f.n3r, we_all.f.n3r, ns_we_all.f.3r, we_all.f.3r),
 				cbind(ns_we_all.f.n4, we_all.f.n4, ns_we_all.f.4, we_all.f.4),
 				cbind(ns_we_all.f.nx, we_all.f.nx, ns_we_all.f.x, we_all.f.x),
 				cbind(ns_we_all.f.ny, we_all.f.ny, ns_we_all.f.y, we_all.f.y))
 			
reg <- c("5'UTR","exon","intron","3'UTR","promoter", "1000bp-upstream", "1000bp-downstream",
		"chr2L", "chr2R","chr3L","chr3R","chr4","chrX","chrY")
write.table(t(c("Feature", "No all non-significant TEs", "No all significant TEs", 
			"No non-significant TEs in feature", "No non-significant TEs in feature", "Pval", "q-val")),
			file=paste(parent.o, "fisher.test_tes_we_features_reed.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
write.table(cbind(reg,counts,pvals, qval), file=paste(parent.o, "fisher.test_tes_we_features_reed.txt", sep=""), 
			col.names=F, row.names=F, quote=F, sep="\t", append=T)

#### Perform the enrichment for TE families
# WE TEs
s <- read.csv(file=paste(parent.o, "significant_we_tes_reed_features_gtf2.txt", sep=""), 
				sep="\t", header=T, stringsAsFactor=F)
ns <- read.csv(file=paste(parent.o, "non-significant_we_tes_reed_features_gtf2.txt", sep=""), 
				sep="\t", header=T, stringsAsFactor=F)

all.r <- all.c <- num.s.f <- num.ns.f <- num.s.nf <- num.ns.nf <- NULL
for (i in 1:length(unique(s[,2]))) {
	all.c <- NULL
	num.s.f <- dim(s[s[,2] == unique(s[,2])[i], ])[1]
	num.ns.f <- dim(ns[ns[,2] == unique(s[,2])[i], ])[1]
	num.s.nf <- dim(s[s[,2] != unique(s[,2])[i], ])[1]
	num.ns.nf <- dim(ns[ns[,2] != unique(s[,2])[i], ])[1]
	all.c <- cbind(num.s.f, num.ns.f, num.s.nf, num.ns.nf)
	all.r <- rbind(all.r, all.c) 
}	
rownames(all.r) <- unique(s[,2])	
colnames(all.r) <- c("significant focus TE", "non-significant focus TE", "significant non-focus TE", 
					"non-significant non-focus TE")

f.fisher <- NULL
for (i in 1:length(all.r[,1])) {
	f.fisher <- c(f.fisher, fisher.test(matrix(all.r[i,], ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value)	
}

all.f <- NULL
qval <- p.adjust(f.fisher, method = "fdr", n = length(f.fisher)) 
all.f <- cbind(all.r, f.fisher, qval)
colnames(all.f)[5:6] <- c("p-val", "q-val")
write.table(all.f, file=paste(parent.o, "fisher.test_we_te_families_reed.txt", sep=""), 
			sep="\t", quote=F, col.names=T, row.names=T)

#### Perform the enrichment for TE families (only families with => 20 TEs per family)
## WE TEs
s <- read.csv(file=paste(parent.o, "significant_we_tes_reed_features_gtf2.txt", sep=""), 
				sep="\t", header=T, stringsAsFactor=F)
ns <- read.csv(file=paste(parent.o, "non-significant_we_tes_reed_features_gtf2.txt", sep=""), 
				sep="\t", header=T, stringsAsFactor=F)

# Get the names of TE families which have => than 20 detected TEs (Te.y). te.n are the families with < 20 copies per family
te.y <- te.n <- NULL 
for (i in 1:length(unique(ns[,2]))) {
	if (dim(ns[ns[,2] == unique(ns[,2])[i], ])[1] >= 20) {
		te.y <- c(te.y, unique(ns[,2])[i])
	} else { te.n <- c(te.n, unique(ns[,2])[i]) }
}		

if (length(te.n) == 0) {
	#Count the numbers of significant and non significant "focus" and "non-focus" TEs
	all.r <- all.c <- num.s.f <- num.ns.f <- num.s.nf <- num.ns.nf <- NULL
	for (i in 1:length(unique(s[,2]))) {
		all.c <- NULL
		num.s.f <- dim(s[s[,2] == unique(s[,2])[i], ])[1]
		num.ns.f <- dim(ns[ns[,2] == unique(s[,2])[i], ])[1]
		num.s.nf <- dim(s[s[,2] != unique(s[,2])[i], ])[1]
		num.ns.nf <- dim(ns[ns[,2] != unique(s[,2])[i], ])[1]
		all.c <- cbind(num.s.f, num.ns.f, num.s.nf, num.ns.nf)
		all.r <- rbind(all.r, all.c) 
	}	
	rownames(all.r) <- unique(s[,2])	
	colnames(all.r) <- c("significant focus TE", "non-significant focus TE", "significant non-focus TE", 
						"non-significant non-focus TE")

	f.fisher <- NULL
	for (i in 1:length(all.r[,1])) {
		f.fisher <- c(f.fisher, fisher.test(matrix(all.r[i,], ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value)	
	}

	qval <- p.adjust(f.fisher, method = "fdr", n = length(f.fisher)) 
	all.f <- cbind(all.r, f.fisher, qval)
	colnames(all.f)[5:6] <- c("p-val", "q-val")
	write.table(all.f, file=paste(parent.o, "fisher.test_we_te_families_reed_gtf2.txt", sep=""), 
				sep="\t", quote=F, col.names=T, row.names=T)
} else {
	s.f <- ns.f <- NULL

	# Get all significant TEs which presence (detected) have more than 20 copies per family 
	for (i in 1:length(te.y)) {
		s.f <- rbind(s.f, s[s[,2] == te.y[i],])
	}

	# Get all non-significant TEs which presence (detected) have more than 20 copies per family 
	for (i in 1:length(te.y)) {
		ns.f <- rbind(ns.f, ns[ns[,2] == te.y[i],])
	}
	
	# Count the numbers of significant and non significant "focus" and "non-focus" TEs
	all.r <- all.c <- num.s.f <- num.ns.f <- num.s.nf <- num.ns.nf <- NULL
	for (i in 1:length(unique(s.f[,2]))) {
		all.c <- NULL
		num.s.f <- dim(s.f[s.f[,2] == unique(s.f[,2])[i], ])[1]
		num.ns.f <- dim(ns.f[ns.f[,2] == unique(s.f[,2])[i], ])[1]
		num.s.nf <- dim(s.f[s.f[,2] != unique(s.f[,2])[i], ])[1]
		num.ns.nf <- dim(ns.f[ns.f[,2] != unique(s.f[,2])[i], ])[1]
		all.c <- cbind(num.s.f, num.ns.f, num.s.nf, num.ns.nf)
		all.r <- rbind(all.r, all.c) 
	}	
	rownames(all.r) <- unique(s.f[,2])	
	colnames(all.r) <- c("significant focus TE", "non-significant focus TE", "significant non-focus TE", 
						"non-significant non-focus TE")

	f.fisher <- NULL
	for (i in 1:length(all.r[,1])) {
		f.fisher <- c(f.fisher, fisher.test(matrix(all.r[i,], ncol = 2, nrow = 2, byrow=T), alternative = "greater")$p.value)	
	}

	qval <- p.adjust(f.fisher, method = "fdr", n = length(f.fisher)) 
	all.f <- cbind(all.r, f.fisher, qval)
	colnames(all.f)[5:6] <- c("p-val", "q-val")
	write.table(all.f, file=paste(parent.o, "fisher.test_we_te_families_more_than_20_reed.txt", sep=""), sep="\t", quote=F, col.names=T, row.names=T)
}