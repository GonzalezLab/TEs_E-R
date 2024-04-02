#!/usr/bin/env Rscript

library(DescTools) 
library(GenomicRanges)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with the TEMP output files (i.e., "/home/")
#args[2]: path to the results folder 

### Read TEMP output files (avoid this step if the script 1.Reed.r was already run before)
parent.f <- args[1]
parent.o <- args[2]
list.files  <- dir(parent.f, pattern ="*summary")

for (d in 1:length(list.files)) {
	setwd(parent.f)
	input <- list.files[d]
	name <- strsplit(input, split="-", fixed=T)[[1]][1]
	output <- paste(name, ".insertion.1p1.f0.10.txt", sep="")

	# Read input	
	x <- read.table(input, header=T, sep="\t", fill=T, row.names=NULL)
	y <- x
	colnames(y) <- colnames(x)[-1]

    # Get only those insertions supported by reads at both sides and with freq > 0.10
	ins.1p1.0.10 <- y[as.numeric(as.character(y[,8])) > 0.10 & y[,6] == "1p1", ] 
	
    # Predicted junctions with bp-resolution (columns 10,11:13 with values different from 0) = Narrow junctions
	ins.1p1.0.10.ref.1 <- NULL
	for (i in 1:length(ins.1p1.0.10[,1])) {
		if (as.numeric(ins.1p1.0.10[i,10]) != 0 && as.numeric(ins.1p1.0.10[i,12]) != 0 && as.numeric(ins.1p1.0.10[i,13]) != 0 && as.numeric(ins.1p1.0.10[i,14]) != 0) {		
			ins.1p1.0.10.ref.1 <- rbind(ins.1p1.0.10.ref.1, ins.1p1.0.10[i, c(1,9,11,4:8,10,12:14)])		
		} 
	}
		
    # Predicted junctions with bp-resolution (columns 10,11:13 with some value = 0) = Broad junctions
	ins.1p1.0.10.ref.2 <- NULL
	for (i in 1:length(ins.1p1.0.10[,1])) {	
		if (as.numeric(ins.1p1.0.10[i,10]) == 0 || as.numeric(ins.1p1.0.10[i,12]) == 0 || as.numeric(ins.1p1.0.10[i,13]) == 0 || as.numeric(ins.1p1.0.10[i,14]) == 0) {	 
			ins.1p1.0.10.ref.2 <- rbind(ins.1p1.0.10.ref.2, ins.1p1.0.10[i, c(1:8,10,12:14)]) 
		}
	}	
		
	# Include the name of the strain for each input file (ex: BL, LS, HS...)	
	strain <- rep(name, length(ins.1p1.0.10.ref.1[,1]))
	ins.1p1.0.10.ref.1 <- cbind(ins.1p1.0.10.ref.1, strain)
	strain <- rep(name, length(ins.1p1.0.10.ref.2[,1]))
	ins.1p1.0.10.ref.2 <- cbind(ins.1p1.0.10.ref.2, strain)	

    # Join both the Narrow and Broad junction and name the final data.frame
	ins.1p1.0.10.ref.1 <- unname(ins.1p1.0.10.ref.1)
	ins.1p1.0.10.ref.2 <- unname(ins.1p1.0.10.ref.2)
	colnames(ins.1p1.0.10.ref.1) <- c("Chr", "Start", "End", "TransposonName", "TransposonDirection", 
										"Class", "VariantSupport", "Frequency", "Junction1Support", 
										"Junction2Support", "5'_Support", "3'_Support", "Strain")
	colnames(ins.1p1.0.10.ref.2) <- c("Chr", "Start", "End", "TransposonName", "TransposonDirection", 
										"Class", "VariantSupport", "Frequency", "Junction1Support", 
										"Junction2Support", "5'_Support", "3'_Support", "Strain")
	ins.1p1.0.10.ref <- rbind(ins.1p1.0.10.ref.1, ins.1p1.0.10.ref.2)
	
    # Remove those TEs that are the same in the same chr and coordinates but annotated twice#1p1	
	data.f0.10.final <- data.f0.10.final.j <- data.f0.10.1.final <- data.f0.10.2.final <- mr.f0.10 <- mr1.f0.10 <- NULL

	# Non-duplicated TEs
	te.j.f0.10 <- ins.1p1.0.10.ref %>% unite(All, c(Chr, Start,  End, TransposonName), sep = "//", remove = TRUE)	

	# Duplicated TEs that are duplicated
	te.d.f0.10 <- te.j.f0.10[duplicated(te.j.f0.10[,1]), ]
	data.f0.10.1.final <- te.j.f0.10[!(te.j.f0.10[,1] %in% te.d.f0.10[,1]), ]

	# Select, among the duplicated TEs, the annotation with more supported reads 
	for (i in 1:length(te.d.f0.10[,1])) {	
	 	ind.f0.10 <- which(te.j.f0.10[,1] == te.d.f0.10[i, 1])
	 	mr.f0.10 <- NULL
	 	for (j in 1:length(ind.f0.10)) {
	 		mr.f0.10 <- rbind(mr.f0.10, te.j.f0.10[ind.f0.10[j], ])
	 	}
	 	mr1.f0.10 <- mr.f0.10[mr.f0.10[,4] == max(mr.f0.10[,4]), ]
	 	data.f0.10.2.final <- rbind(data.f0.10.2.final, mr1.f0.10[1, ])
	}	 

	# Join the non-duplicated TEs and those among the duplicated with more supponting reads
	data.f0.10.final.j <- rbind(data.f0.10.1.final, data.f0.10.2.final)	 
	data.f0.10.final.j <- na.omit(data.f0.10.final.j) 
	data.f0.10.final <- cbind(matrix(unlist(strsplit(data.f0.10.final.j$All, split="//", fixed=T)), 
								ncol=4, byrow=T, dimnames=list(NULL, c("Chr", "Start", "End", "TransposonName"))), 
								data.f0.10.final.j[, c(2:10)])
	data.f0.10.final[, 2] <- as.numeric(as.character(data.f0.10.final[,2]))
	data.f0.10.final[, 3] <- as.numeric(as.character(data.f0.10.final[,3]))
	
	write.table(data.f0.10.final, file=paste(parent.o, output, sep=""), quote=F, col.names=T, row.names=F, sep="\t")	 	
}

### Compare each CONTROL and SELECTION pool to get the frequencies of the intersected intervals and perform Fisher's exact test and Bonferroni correction
list.files.temp  <- dir(parent.o, pattern ="*insertion.1p1.f0.10.txt")

for (f in 1:length(list.files.temp)) {
	temp.f <- NULL
	input.temp <- paste(parent.o, list.files.temp[f], sep="")
	temp <- read.table(file=paste(parent.o, list.files.temp[f], sep=""), header=T)
	temp.f <- cbind(temp, matrix("NA", ncol=1, nrow=length(temp[,1])))
	temp.f[, 14] <- as.numeric(temp.f[, 14])
	name <- strsplit(list.files.temp[f], split=".", fixed=T)[[1]][1]
	
	for (l in 1:length(temp[,1])) {
		ref <- round((as.numeric(temp[l,7]) - (as.numeric(temp[l,8])*as.numeric(temp[l,7])))/as.numeric(temp[l,8]))	
		temp.f[l, 14] <- round((as.numeric(temp[l,7]) - (as.numeric(temp[l,8])*as.numeric(temp[l,7])))/as.numeric(temp[l,8]))	
	}
	colnames(temp.f)[14]<- "Reads_ref"
	write.table(temp.f, file=paste(parent.o, name, ".insertion.1p1.f0.10.final_ref.txt", sep=""), quote=F, 
				sep="\t", row.names=F, col.names=T)
}

list.files.bl.1  <- dir(parent.o, pattern ="BL.*insertion.1p1.f0.10.final_ref.txt")

# Append the intervals for all BL pools
all.bs.1 <- NULL
for (i in 1:length(list.files.bl.1)) {
	bs.1 <- read.table(file=paste(parent.o,list.files.bl.1[i], sep=""), header=T, sep="\t") 
	all.bs.1 <- rbind(all.bs.1, bs.1)
}

# Get new intervals from the appended BL pools to remove overlapping.
all.bs.1u <- all.bs.1 %>% unite(Chr_TransposonName, c(Chr, TransposonName), sep = "&", remove = TRUE)
all.bs.1r <- GRanges(all.bs.1u$Chr_TransposonName, IRanges(all.bs.1u$Start, all.bs.1u$End))
all.bs.1r.red <- reduce(all.bs.1r)
all.1 <- as.data.frame(all.bs.1r.red)
all.1.s <- all.1[,1:3]
colnames(all.1.s) <- c("Chr_TransposonName", "Start", "End")

# Compare the new intervals that do not overlap with each BL pool and get the interval among the BL pools with maximum reads 
dat <- r.bl <- NULL
for (i in 1:length(all.1.s[,1])) { 
	dat <- NULL
	for (j in 1:length(list.files.bl.1)) {
		x <- read.table(file=paste(parent.o,list.files.bl.1[j], sep=""), header=T, sep="\t") 
		x.u <- x %>% unite(Chr_TransposonName, c(Chr, TransposonName), sep = "&", remove = TRUE)
		for (k in 1:length(x.u[,1])) { 
			if (all.1.s[i,1] == x.u[k,1]) {
				if (c(all.1.s[i,2], all.1.s[i,3]) %overlaps% c(x.u[k,2], x.u[k,3])) {
					dat <- rbind(dat, x.u[k, ]) 
				}	
			}				
		}
	}
	if (!is.null(dat) && dim(dat)[1] > 1) { 
		m <- max(rowSums(dat[,c(6,13)]))
		su <- rowSums(dat[,c(6,13)])
		dat.max <- dat[names(su[su == m]), ] 
		if ( dim(dat.max)[1] > 1 ) { 
		 	if (dim(dat.max[dat.max[,7] == max(dat.max[,7]), ])[1] == 1) { 
		 		r.bl <- rbind(r.bl, dat.max[dat.max[,7] == max(dat.max[,7]), ])
		 	} else {
		 		r.bl <- rbind(r.bl,  dat.max[dat.max[,7] == max(dat.max[,7]), ][1,])
		 	}
		} else {
			r.bl <- rbind(r.bl, dat.max)
		}
	}
	if ( !is.null(dat) && dim(dat)[1] == 1) {
		r.bl <- rbind(r.bl, dat)
	}	 	
}	

r.bl.final <- cbind(matrix(unlist(strsplit(as.character(r.bl$Chr_TransposonName), split="&", fixed=T)), 
					ncol=2, byrow=T, dimnames=list(NULL, c("Chr", "TransposonName"))),r.bl[,2:13]) 

write.table(all.1.s, file=paste(parent.o, "concatenated.all.bl.1p1.f0.10.txt", sep=""), quote=F, sep="\t", 
			row.names=F, col.names=T) 
write.table(r.bl.final, file=paste(parent.o, "concatenated.all.bl.1p1.f0.10.final.txt", sep=""), quote=F, 
			sep="\t", row.names=F, col.names=T) 

# Compare BL with each SELECTION pool to get the frequencies and perform Fisher's exact test and Bonferroni correction
input.c <- read.table(file=paste(parent.o, "concatenated.all.bl.1p1.f0.10.final.txt", sep=""), sep="\t", header=T)

# HS pool
for (f in 1:2) {	
	input.r <- read.table(file=paste(parent.o, "HS",f,".insertion.1p1.f0.10.final_ref.txt", sep=""), sep="\t", header=T)
	data1 <- NULL
	for (i in 1:length(input.c[,1])) {
		for (j in 1:length(input.r[,1])) {
			if (input.c[i,1] == input.r[j,1] && input.c[i,2] == input.r[j,4]) {
				if (c(input.c[i,3],input.c[i,4]) %overlaps% c(input.r[j,2], input.r[j,3])) {
					pval <- fisher.test(matrix(c(as.numeric(input.c[i,7]), as.numeric(input.c[i,14]), 
										as.numeric(input.r[j,7]), as.numeric(input.r[j,14])), nrow=2, ncol=2, byrow=T))$p.val
					data1 <- rbind(data1, c(as.character(input.c[i,1]), as.character(input.c[i,2]), 
									as.numeric(input.c[i,3]), as.numeric(input.c[i,4]), as.numeric(input.r[j,2]), as.numeric(input.r[j,3]), as.numeric(input.c[i,8]), as.numeric(input.r[j,8]), scientific(as.numeric(pval), digits=2)))
				}
			}		
		}
	}
	colnames(data1) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval")
	write.table(data1, file=paste(parent.o, "Freqs_BL-HS", f,"_f0.10_Pval_reads.txt", sep=""), quote=F, sep="\t", 
				row.names=F, col.names=T)

	data2 <- data1
	data2 <- cbind(data2, rep("-", length(data2[,1])) )
	for (d in 1:length(data2[,1])) {
		if ( as.numeric(data2[d,9]) < (0.05/length(data2[,1])) ) { 
			data2[d,10] <- "*"
		}
		if ( as.numeric(data2[d,9]) < (0.01/length(data2[,1])) ) { 
			data2[d,10] <- "**"
		}
		colnames(data2) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval", 
							"Pval_BF")
		write.table(data2, file=paste(parent.o, "Freqs_BL-HS", f,"_f0.10_Pval_BF_reads.txt", sep=""), quote=F, 
					sep="\t", row.names=F, col.names=T)
	}
}

# LS pool
for (f in 1:2) {	
	input.r <- read.table(file=paste(parent.o, "LS",f,".insertion.1p1.f0.10.final_ref.txt", sep=""), sep="\t", header=T)
	data1 <- NULL
	for (i in 1:length(input.c[,1])) {
		for (j in 1:length(input.r[,1])) {
			if (input.c[i,1] == input.r[j,1] && input.c[i,2] == input.r[j,4]) {
				if (c(input.c[i,3],input.c[i,4]) %overlaps% c(input.r[j,2], input.r[j,3])) {
					pval <- fisher.test(matrix(c(as.numeric(input.c[i,7]), as.numeric(input.c[i,14]), 
										as.numeric(input.r[j,7]), as.numeric(input.r[j,14])), nrow=2, ncol=2, byrow=T))$p.val
					data1 <- rbind(data1, c(as.character(input.c[i,1]), as.character(input.c[i,2]), as.numeric(input.c[i,3]), 
									as.numeric(input.c[i,4]), as.numeric(input.r[j,2]), as.numeric(input.r[j,3]), as.numeric(input.c[i,8]),
									as.numeric(input.r[j,8]), scientific(as.numeric(pval), digits=2)))
				}
			}		
		}
	}
	colnames(data1) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval")
	write.table(data1, file=paste(parent.o, "Freqs_BL-LS", f,"_f0.10_Pval_reads.txt", sep=""), quote=F, sep="\t",
									row.names=F, col.names=T)

	data2 <- data1
	data2 <- cbind(data2, rep("-", length(data2[,1])) )
	for (d in 1:length(data2[,1])) {
		if ( as.numeric(data2[d,9]) < (0.05/length(data2[,1])) ) { 
			data2[d,10] <- "*"
		}
		if ( as.numeric(data2[d,9]) < (0.01/length(data2[,1])) ) { 
			data2[d,10] <- "**"
		}
		colnames(data2) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval", 
							"Pval_BF")
		write.table(data2, file=paste(parent.o, "Freqs_BL-LS", f,"_f0.10_Pval_BF_reads.txt", sep=""), quote=F, 
					sep="\t", row.names=F, col.names=T)
	}
}

# HF pool
for (f in 1:2) {	
	input.r <- read.table(file=paste(parent.o, "HF",f,".insertion.1p1.f0.10.final_ref.txt", sep=""), sep="\t", header=T)
	data1 <- NULL
	for (i in 1:length(input.c[,1])) {
		for (j in 1:length(input.r[,1])) {
			if (input.c[i,1] == input.r[j,1] && input.c[i,2] == input.r[j,4]) {
				if (c(input.c[i,3],input.c[i,4]) %overlaps% c(input.r[j,2], input.r[j,3])) {
					pval <- fisher.test(matrix(c(as.numeric(input.c[i,7]), as.numeric(input.c[i,14]), 
										as.numeric(input.r[j,7]), as.numeric(input.r[j,14])), nrow=2, ncol=2, byrow=T))$p.val
					data1 <- rbind(data1, c(as.character(input.c[i,1]), as.character(input.c[i,2]), 
									as.numeric(input.c[i,3]), as.numeric(input.c[i,4]), as.numeric(input.r[j,2]), 
									as.numeric(input.r[j,3]), as.numeric(input.c[i,8]), as.numeric(input.r[j,8]), 
									scientific(as.numeric(pval), digits=2)))
				}
			}		
		}
	}
	colnames(data1) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval")
	write.table(data1, file=paste(parent.o, "Freqs_BL-HF", f,"_f0.10_Pval_reads.txt", sep=""), 
									quote=F, sep="\t", row.names=F, col.names=T)

	data2 <- data1
	data2 <- cbind(data2, rep("-", length(data2[,1])) )
	for (d in 1:length(data2[,1])) {
		if ( as.numeric(data2[d,9]) < (0.05/length(data2[,1])) ) { 
			data2[d,10] <- "*"
		}
		if ( as.numeric(data2[d,9]) < (0.01/length(data2[,1])) ) {
			data2[d,10] <- "**"
		}
		colnames(data2) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval", 
							"Pval_BF")
		write.table(data2, file=paste(parent.o, "Freqs_BL-HF", f,"_f0.10_Pval_BF_reads.txt", sep=""), quote=F,
					 sep="\t", row.names=F, col.names=T)
	}
}