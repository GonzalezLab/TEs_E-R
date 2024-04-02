#!/usr/bin/env Rscript

library(DescTools) 
library(GenomicRanges)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with containing the TEMP input files obtained from running the script Griffin.r (i.e., "/home/")
#args[2]: path to the results folder

parent.f <- args[1]
parent.o <- args[2]

### Merge all the TE intervals from all the pools
list.files.1  <- list.files(parent.f, pattern = "*_new_intervals_1p1.f0.10_reads.txt")

int.ref.1 <-  NULL
for (f in 1:length(list.files.1)) {
	input.1 <- paste(parent.f, list.files.1[f], sep="") 
	data.1 <- read.table(file=input.1, header=T)
	data.1 <- cbind(data.1, rep(strsplit(list.files.1[f], "_", fixed=T)[[1]][1], length(data.1[,1])))
	colnames(data.1)[8] <- "pool"
	int.ref.1 <- rbind(int.ref.1, data.1)
}
	
write.table(int.ref.1, file=paste(parent.o, "griffin_intervals_concatenated_1p1.f0.10.txt", sep=""), 
			quote=F, col.names=T, row.names=F, sep="\t")

### Get the intervals for the non-significant TEs
conc <- read.table(file=paste(parent.o, "griffin_intervals_concatenated_1p1.f0.10.txt", sep=""), 
					sep="\t", header=T)

# Merge overlapping intervals in each control/selection pool
conc.n <- NULL
for (i in 1:length(unique(conc[,8]))) {
	p <- conc[conc[,8] == unique(conc[,8])[i], ]
	strain.name <- unique(conc[,8])[i]
	p.p <- p %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
	p.ranges <- GRanges(p.p$Chr_TE_ID, IRanges(p.p$Start, p.p$Stop))
	p.ranges.red <- reduce(p.ranges)
	p.new.int <- as.data.frame(p.ranges.red)
	p.new <- cbind(p.new.int[, 1:3], rep(p[1,8], length(p.new.int[,1])))
	colnames(p.new) <- c("Chr_TE_ID", "Start", "Stop", "Comparison")
	conc.n <- rbind(conc.n, p.new)
}

# Merge overlapping intervals from all control/selection pool 	
conc.f <- GRanges(conc.n$Chr_TE_ID, IRanges(conc.n$Start, conc.n$Stop))
conc.f.red <- reduce(conc.f)
nred <- as.data.frame(conc.f.red)
nred <- nred[ ,1:3]
colnames(nred) <- c("Chr_TE_ID", "Start", "Stop")

int.new.1 <- cbind(nred, rep(0, length(nred[,1])))
colnames(int.new.1)[4] <- "Freq"

for (f in 1:length(list.files.1)) {
	input.1 <- paste(parent.f, list.files.1[f], sep="")
	data.1.1 <- read.table(file=input.1, header=T)
	if (!is.null(data.1.1) && dim(data.1.1) != 0) {
		data.1 <- NULL
		name <- strsplit(list.files.1[f], split="_", fixed=T)[[1]][1]
		id1 <- id2 <- id3 <- id4 <- NULL
		id1 <- rep("ID", length(data.1.1[,1]))
		id2 <- seq(1, length(data.1.1[,1]), 1)
		id3 <- rep(name, length(data.1.1[,1]))
		id4 <- paste(id3,"_", id1, id2, sep="")
		data.1 <- cbind(data.1.1, id4)
		colnames(data.1)[7] <- "Strain_ID"
		
		write.table(data.1, file=paste(parent.o, strsplit(list.files.1[f], ".", fixed=T)[[1]][1], "_ID.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")
		data.1 <- data.1 %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
		for (i in 1:length(int.new.1[,1])) {
			for (j in 1:length(data.1[,1])) {
				if ( int.new.1[i,1] == data.1[j,1] ) {
					if (c(int.new.1[i,2], int.new.1[i,3]) %overlaps% c(data.1[j,2], data.1[j,3]) == TRUE) {	
						int.new.1[i,4] <- int.new.1[i,4] + 1 
				 		if (is.null(int.new.1[i,5]) || is.na(int.new.1[i,5])) {
				 			int.new.1[i,5] <- as.character(data.1[j,7])
				 		} else {
				 		int.new.1[i,5] <- paste(int.new.1[i,5], as.character(data.1[j,7]), sep=";", collapse=NULL)
				 		}
				 	}
				}
			}	 	
		}
	}
}

colnames(int.new.1)[5] <- "Pool"
write.table(int.new.1, file=paste(parent.o, "summary_intervals_and_TE_frequency_1p1.f0.10.txt", 
sep=""), quote=F, col.names=T, row.names=F, sep="\t")

x <- read.table(file=paste(parent.o, "summary_intervals_and_TE_frequency_1p1.f0.10.txt", sep=""), 
				header=T, sep="\t", stringsAsFactors=F)
sig <- x[x[,4] == 1 & grepl("DR", x[,5], fixed=T), ]
con <- x[x[,4] == 1 & grepl("CR", x[,5], fixed=T), ]

list.files.s  <- list.files(parent.o, pattern = "DR.*_new_intervals_1p1_ID.txt")
list.files.c  <- list.files(parent.o, pattern = "CR.*_new_intervals_1p1_ID.txt")

sig.f <- NULL
con.f <- NULL		

### Get the frequency of those TEs that are only in one DR pool
for (i in 1:length(list.files.s)) {
	for (j in 1:length(sig[,1])) {
		y <- read.table(file=paste(parent.o, list.files.s[i], sep=""), header=T,sep="\t", 
							stringsAsFactors=F) 
		for (z in 1:length(y[,1])) {
			if (sig[j ,5] ==  y[z, 8]) {
				sig.f <- rbind(sig.f, y[z, ])
			}
		}
	}		
}
colnames(sig.f)[8] <- "ID"	
write.table(sig.f, file=paste(parent.o, "summary_TEs_single_DR_frequency_1p1.f0.10.txt", sep=""), 
							quote=F, col.names=T, row.names=F, sep="\t")
							
### Get the frequency of the TEs that are only in one CR pool	
for (i in 1:length(list.files.c)) {
	for (j in 1:length(con[,1])) {
		y <- read.table(file=paste(parent.o, list.files.c[i], sep=""), header=T,sep="\t", 
							stringsAsFactors=F) 
		for (z in 1:length(y[,1])) {
			if (con[j ,5] ==  y[z, 8]) {
				con.f <- rbind(con.f, y[z, ])
			}
		}
	}		
}
colnames(con.f)[8] <- "ID"	
write.table(con.f, file=paste(parent.o, "summary_TEs_single_CR_frequency_1p1.f0.10.txt", sep=""), 
				quote=F, col.names=T, row.names=F, sep="\t")

# Test if frequencies in DR are significantly different from frequencies in CR (Wilcoxon test)
i.dr <- "summary_TEs_single_DR_frequency_1p1.f0.10.txt"
i.cr <- "summary_TEs_single_CR_frequency_1p1.f0.10.txt"
x <- read.table(file=paste(parent.o, i.dr, sep=""), header=T, sep="\t", stringsAsFactors=F)
y <- read.table(file=paste(parent.o, i.cr, sep=""), header=T, sep="\t", stringsAsFactors=F)

wilcox.test(x[,5], y[,5], alternative = "greater"))
