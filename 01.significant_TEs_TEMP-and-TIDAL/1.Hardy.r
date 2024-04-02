#!/usr/bin/env Rscript

library(DescTools) 
library(GenomicRanges)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with the diferent subfolders (pools) containing the TEMP output files (i.e., "/home/")
#args[2]: path to the results folder
#args[3]: file with the TEs dictionary FlyBase-Repbase (with full path)
#args[4]: path to the folder with the diferent subfolders (pools) containing the TIDAL output files (i.e., "/home/")

### Read TEMP output files
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
		if (as.numeric(ins.1p1.0.10[i,10]) != 0 && as.numeric(ins.1p1.0.10[i,12]) != 0 && as.numeric(ins.1p1.0.10[i,13])!= 0 && as.numeric(ins.1p1.0.10[i,14]) != 0) {		
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
		
	# Include the name of the strain for each input file (ex: CA, CB...)	
	strain <- rep(name, length(ins.1p1.0.10.ref.1[,1]))
	ins.1p1.0.10.ref.1 <- cbind(ins.1p1.0.10.ref.1, strain)
	strain <- rep(name, length(ins.1p1.0.10.ref.2[,1]))
	ins.1p1.0.10.ref.2 <- cbind(ins.1p1.0.10.ref.2, strain)	

	# Join both the Narrow and Broad junction and name the final data.frame
	ins.1p1.0.10.ref.1 <- unname(ins.1p1.0.10.ref.1)
	ins.1p1.0.10.ref.2 <- unname(ins.1p1.0.10.ref.2)
	colnames(ins.1p1.0.10.ref.1) <- c("Chr", "Start", "End", "TransposonName", "TransposonDirection", "Class", 
										"VariantSupport", "Frequency", "Junction1Support", "Junction2Support", 
										"5'_Support", "3'_Support", "Strain")
	colnames(ins.1p1.0.10.ref.2) <- c("Chr", "Start", "End", "TransposonName", "TransposonDirection", "Class", 
										"VariantSupport", "Frequency", "Junction1Support", "Junction2Support", 
										"5'_Support", "3'_Support", "Strain")
	ins.1p1.0.10.ref <- rbind(ins.1p1.0.10.ref.1, ins.1p1.0.10.ref.2)
	
	# Remove those TEs that are the same in the same chr and coordinates but annotated twice
	data.f0.10.final <- data.f0.10.final.j <- data.f0.10.1.final <- data.f0.10.2.final <- mr.f0.10 <- mr1.f0.10 <- NULL

	# Non-duplicated TEs	
	te.j.f0.10 <- ins.1p1.0.10.ref %>% unite(All, c(Chr, Start,  End, TransposonName), sep = "//", remove = TRUE)	

	# Duplicated TEs	
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

### Change nomenclature in TEMP (REPBASE nomenclature for TEs) files for uniformity with TIDAL (FlyBase nomenclature)
list.files  <- dir(parent.o, pattern ="*insertion.1p1.*txt")
nom.file <- args[3]
nom <- read.table(nom.file, header=T, sep="\t")

for (f in 1:length(list.files)) {
 		x <- read.table(file=paste(parent.o,list.files[f], sep=""), header=T, sep="\t", fill=T, row.names=NULL)
 		colnames(x) <- c("Chr", "Start", "End", "TransposonName", "TransposonDirection", "Class", "VariantSupport", 
 						"Frequency", "Junction1Support", "Junction2Support", "p5_Support", "p3_Support", "Strain")
 		y <- x
		y[,4] <- as.character(y[,4])
		for (l in 1:length(y[,1])) {
			if (as.character(y[l,4]) %in% as.character(nom[,2])) {
				y[l,4] <- as.character(unique(nom[as.character(nom[ ,2]) == as.character(y[l,4]), 1])[1]) #Assign the first name of FB that match the with the REPBASE name in TEMP file
			} else {	y[l,4] <- y[l,4] }
 		}
 		nam <-strsplit(list.files[f], split=".txt", fixed=T)[[1]][1]
 
		y.join <- y %>% unite(All, c(Chr, TransposonName), sep = "//", remove = TRUE)
		p.ranges <- GRanges(y.join$All, IRanges(y.join$Start, y.join$End))
		p.ranges.red <- reduce(p.ranges)
		p.int <- as.data.frame(p.ranges.red)

		p <- p.int[ ,1:3]
		colnames(p) <- c("All", "Start", "End")
		ref <- p 
		new.y <- y %>% unite(All, c(Chr, TransposonName), sep = "//", remove = TRUE)

		a <- b <- c <- cc <- NULL
		for (i in 1:length(ref[,1])) {
			a <- NULL
			for (j in 1:length(new.y[,1])) {				
				if (ref[i,1] == new.y[j,1]) {
					if (c(ref[i,2],ref[i,3]) %overlaps% c(new.y[j,2],new.y[j,3])) {
						a <- rbind(a, new.y[j, ])
					} 
				}
			}
			b <- a[as.numeric(a[,4]) == max(as.numeric(a[,4])), ]
			c <- rbind(c, b[1, ])			
		}

		cc <- matrix(unlist(strsplit(c$All, split="//", fixed=T)), ncol=2, byrow=T, 
					dimnames=list(NULL, c("Chr", "TransposonName")))
		m.final <- cbind(cc[ ,1], c[ ,2], c[ ,3], cc[ ,2], c[ ,4:12])
		colnames(m.final)[1:4] <- c("Chr", "Start", "End", "TransposonName")
		write.table(m.final, paste(parent.o, nam, ".final.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")
}	
  
### Read TIDAL output files
tidal.folder <- args[4]

list.files.temp  <- dir(parent.o, sep=""), pattern ="*insertion.1p1.f0.10.final.txt")
list.files.tidal  <- dir(tidal.folder, pattern ="*Inserts_Annotated.txt")

### Compare TEMP and TIDAL TE insertion intervals
for (f in 1:length(list.files.temp)) {
	input.temp <- paste(parent.o, list.files.temp[f], sep="")
	input.tidal <- paste(tidal.folder, list.files.tidal[f], sep="")
	i.temp <- strsplit(list.files.temp[f], split=".", fixed=T)[[1]][1]
	i.tidal <- strsplit(list.files.tidal[f], split="_", fixed=T)[[1]][1]
	if ( i.temp == i.tidal ) {
		temp <- read.table(file=paste(parent.o, list.files.temp[f], sep=""), header=T)
		tidal <- read.csv(file=paste(tidal.folder, list.files.tidal[f], sep=""), sep="\t", header=T)
		name <- strsplit(input.temp, ".", fixed=T)[[1]][1]
		
		data1.p <- data1.j <- data1.u <- data1 <- NULL
		for (i in 1:length(temp[,1])) {
			for (j in 1:length(tidal[,1])) {
				if (temp[i,1] == tidal[j,2] & temp[i,4] == tidal[j,5]) {
					if (c(temp[i,2],temp[i,3]) %overlaps% c(tidal[j,3], tidal[j,4])) {
						min.range <- min(c(temp[i,2], temp[i,3], tidal[j,3], tidal[j,4]))
						max.range <- max(c(temp[i,2], temp[i,3], tidal[j,3], tidal[j,4]))
						data1.p <- rbind(data1.p, c(as.character(temp[i,1]), min.range, max.range, as.character(temp[i,4]), 
										temp[i,7], temp[i,8]))
						}
				}		
			}
		}
		colnames(data1.p) <- c("Chr", "Start", "Stop", "TE_ID", "Reads", "Freq")
		
		# Remove redundancy
		data1.j <- as.data.frame(data1.p) %>% unite(All, c(Chr, Start,  Stop, TE_ID, Reads, Freq), sep = "&", remove = TRUE)
		data1.u <- unique(data1.j)
		data1 <- matrix(unlist(strsplit(data1.u$All, split="&", fixed=T)), ncol=6, byrow=T, dimnames=list(NULL, c("Chr", "Start", "Stop", "TE_ID", "Reads", "Freq")))

		# Get the number of reads for the reference alleles
		data1n <- cbind(data1n, rep(NA, length(data1n[,1])))
		for (l in 1:length(data1[,1])) {
			data1n[l ,7] <-  round((as.numeric(data1[l,5]) - (as.numeric(data1[l,6])*as.numeric(data1[l,5])))/as.numeric(data1[l,6]))
		}
		data1n2 <- NULL
		data1n2 <- cbind(data1n[, 1:4],data1n[, 6], round(as.numeric(data1n[, 5])), round(as.numeric(data1n[, 7])))
		colnames(data1n2)[5:7] <- c("Freq", "Reads_alt", "Reads_ref")
		
		data1n2 <- as.data.frame(data1n2)
		data1n2[,1] <- as.character(data1n2[,1]) 
		data1n2[,2] <- as.numeric(as.character(data1n2[,2]))
		data1n2[,3] <- as.numeric(as.character(data1n2[,3]))
		data1n2[,4] <- as.character(data1n2[,4])
		data1n2[,5] <- as.numeric(as.character(data1n2[,5]))
		data1n2[,6] <- as.numeric(as.character(data1n2[,6]))
		data1n2[,7] <- as.numeric(as.character(data1n2[,7]))

		dl <- data1n2[duplicated(data1n2[,2]) & duplicated(data1n2[,1]) | duplicated(data1n2[,3]) & duplicated(data1n2[,1]), ]
		if (dim(dl)[1] != 0 ) {
			mdu <- mjoin.1 <- mrall <- NULL
			for (i in 1:length(dl[,1])) {
				ind <- which(data1n2[ ,1] == dl[i ,1] & data1n2[ ,2] == dl[i ,2] & data1n2[ ,4] == dl[i ,4] | 
							data1n2[ ,1] == dl[i ,1] & data1n2[ ,3] == dl[i ,3] & data1n2[ ,4] == dl[i ,4])
				mr <- NULL
				for (j in 1:length(ind)) {
	 				mr <- rbind(mr, data1n2[ind[j], ])
				}
				mrall <- rbind(mrall, mr)
				mr1 <- mr[as.numeric(as.character(mr[,6])) + as.numeric(as.character(mr[,7])) == max(as.numeric(as.character(mr[,6])) + 
						as.numeric(as.character(mr[,7]))) , ]
				mdu <- rbind(mdu, mr1[1, ])
			}
			colnames(mdu) <- colnames(data1n2)
		mndu <- anti_join(as.data.frame(data1n2, stringsAsFactors=FALSE), as.data.frame(mrall,stringsAsFactors=FALSE))
		m.join.1 <- rbind(mndu, mdu)
		} else {
			m.join.1 <- data1n2	
		}
		
		dl2 <- NULL
		for (i in 1:length(m.join.1[,1])) {
			if( i < length(m.join.1[,1])) {
			for (j in (i+1):length(m.join.1[,1])) {
					if ( m.join.1[i,1] == m.join.1[j,1] & m.join.1[i,4] == m.join.1[j,4] ) {
						if ( c(m.join.1[i,2], m.join.1[i,3]) %overlaps% c(m.join.1[j,2], m.join.1[j,3]) ) {
					 		maxr <- max(as.numeric(as.character(m.join.1[i,6])) + as.numeric(as.character(m.join.1[i,7])), 
					 					as.numeric(as.character(m.join.1[j,6])) + as.numeric(as.character(m.join.1[j,7])))
					 		if (as.numeric(as.character(m.join.1[i,6])) + as.numeric(as.character(m.join.1[i,7])) == maxr) {
					 			dl2 <- rbind(dl2, m.join.1[j, ])
					 		} else {
					 			dl2 <- rbind(dl2, m.join.1[i, ])
					 		}
					 	}
					}
				}
			}	
		}			 		
		if ( !is.null(dim(dl2))[1] ) {
			mndu2.1 <- NULL		
			mndu2.1 <- anti_join(as.data.frame(m.join.1, stringsAsFactors=FALSE), as.data.frame(dl2, stringsAsFactors=FALSE))
			mndu2.1 <- apply(mndu2.1, 2, as.character)	
		} else {
			m.join.1 <- apply(m.join.1,2, as.character)
			mndu2.1 <- m.join.1
		}	
				
		write.table(mndu2.1, file=paste(parent.o, i.temp, "_new_intervals_1p1.f0.10_reads.txt", sep=""), 
					quote=F, sep="\t", row.names=F, col.names=T)
	}
}

### Compare selection pool to each control pool
list.files.sel  <- dir(parent.o, pattern ="S.*new_intervals_1p1.f0.10_reads.txt") 
list.files.con  <- dir(parent.o, pattern ="C.*new_intervals_1p1.f0.10_reads.txt")

for (s in 1:length(list.files.sel)) {
	input.s <- read.table(file=paste(parent.o, list.files.sel[s], sep=""), sep="\t", header=T)
	input.s.n <- input.s %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
    name.s <- strsplit(list.files.sel[s], split="_", fixed=T)[[1]][1]			
	
	# Controls
	for (f in 1:length(list.files.con)) {										
		input.c <- read.table(file=paste(parent.o, list.files.con[f], sep=""), sep="\t", header=T)
		input.c.n <- input.c %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
		data1 <- NULL
		name.c <- strsplit(list.files.con[f], split="_", fixed=T)[[1]][1]
		for (i in 1:length(input.s.n[,1])) {
			for (j in 1:length(input.c.n[,1])) {
				if (input.s.n[i,1] == input.c.n[j,1]) {
					if (c(input.s.n[i,2],input.s.n[i,3]) %overlaps% c(input.c.n[j,2], input.c.n[j,3])) {
						pval <- fisher.test(matrix(c(as.numeric(input.s.n[i,5]), as.numeric(input.s.n[i,6]), 
											as.numeric(input.c.n[j,5]), as.numeric(input.c.n[j,6])), nrow=2, ncol=2, byrow=T))$p.val
						data1 <- rbind(data1, c(strsplit(input.s.n[i,1], split="&&")[[1]][1], strsplit(input.s.n[i,1], split="&&")[[1]][2], 
										input.s.n[i,2], input.s.n[i,3], input.c.n[j,2], input.c.n[j,3], input.s.n[i,4], input.c.n[j,4], 
										scientific(as.numeric(pval), digits=2)))
					}
				}		
			}
		}
		
		# P-values
		colnames(data1) <- c("Chr", "TE_ID", "Start_SEL", "Stop_SEL", "Start_CTR", "Stop_CTR", "Freq_SEL", "Freq_CTR", "Pval")
		write.table(data1, file=paste(parent.o, "Freqs_", name.s, "-", name.c, "_f0.10.Pval_reads.txt", sep=""), 
					quote=F, sep="\t", row.names=F, col.names=T)

		data2 <- data1
		data2 <- cbind(data2, rep("-", length(data2[,1])) )
		for (d in 1:length(data2[,1])) {
			if ( as.numeric(data2[d,9]) < (0.05/length(data2[,1])) ) {
				data2[d,10] <- "*"
			}
			if (as.numeric(data2[d,9]) < (0.01/length(data2[,1])) ) {
				data2[d,10] <- "**"
			}
		}
		
		# Bonferroni Correction
		colnames(data2) <- c("Chr", "TE_ID", "Start_SEL", "Stop_SEL", "Start_CTR", "Stop_CTR", "Freq_SEL", "Freq_CTR", 
							"Pval", "Pval_BF")
		write.table(data2, file=paste(parent.o, "Freqs_", name.s, "-", name.c, "_f0.10.Pval_BF_reads.txt", sep=""), 
					quote=F, sep="\t", row.names=F, col.names=T)
	}
}

### Compare each selection pool to all controls
# Selection pool A
list.files.sel  <- dir(parent.o, pattern ="Freqs_SA-C.*_f0.10.Pval_BF_reads.txt")

input.sa <- read.table(file=paste(parent.o, list.files.sel[1], sep=""), sep="\t", header=T)
input.sb <- read.table(file=paste(parent.o, list.files.sel[2], sep=""), sep="\t", header=T)
input.sc <- read.table(file=paste(parent.o, list.files.sel[3], sep=""), sep="\t", header=T)

sall <- r.sel <- sa <- NULL
sall <- rbind(input.sa, input.sb, input.sc)
sall.u <- sall %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
p.ranges <- GRanges(sall.u$Chr_TE_ID, IRanges(sall.u$Start_SEL, sall.u$Stop_SEL))
p.ranges.red <- reduce(p.ranges)
p.new.int <- as.data.frame(p.ranges.red)
r.sel <- cbind(unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 1)), 
				unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 2)), p.new.int[,2], p.new.int[,3])  
r.sel1 <- as.data.frame(r.sel, stringsAsFactors=FALSE) 
r.sel1[,3] <- as.numeric(r.sel1[,3])
r.sel1[,4] <- as.numeric(r.sel1[,4])
colnames(r.sel1) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new")

sa <- cbind(r.sel1, as.data.frame(matrix("NA", nrow=length(p.new.int[,1]), ncol=24), stringsAsFactors=FALSE))  

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sa[,1])) {
		if (r.sel1[i,1] == input.sa[j,1] && r.sel1[i,2] == input.sa[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sa[j,5], input.sa[j,6])) {
				sa[i, c(5:12)] <- input.sa[j, c(3:10)]
			}
	    }
	}	
}				

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sb[,1])) {
		if (r.sel1[i,1] == input.sb[j,1] && r.sel1[i,2] == input.sb[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sb[j,5], input.sb[j,6])) {
				sa[i, c(13:20)] <- input.sb[j, c(3:10)]
			}
	    }
	}	
}

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sc[,1])) {
		if (r.sel1[i,1] == input.sc[j,1] && r.sel1[i,2] == input.sc[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sc[j,5], input.sc[j,6])) {
				sa[i, c(21:28)] <- input.sc[j, c(3:10)]
			}
	    }
	}	
}

sa[, 12] <- gsub("1", "-", sa[, 12])
sa[, 12] <- gsub("2", "*", sa[, 12])
sa[, 12] <- gsub("3", "**", sa[, 12])
sa[, 20] <- gsub("1", "-", sa[, 20])
sa[, 20] <- gsub("2", "*", sa[, 20])
sa[, 20] <- gsub("3", "**", sa[, 20])
sa[, 28] <- gsub("1", "-", sa[, 28])
sa[, 28] <- gsub("2", "*", sa[, 28])
sa[, 28] <- gsub("3", "**", sa[, 28])

colnames(sa)[5:28] <- c("Start_SELA","Stop_SELA", "Start_CTRA","Stop_CTRA", "Freq_SELA", "Freq_CTRA", "Pval", "Pval_BF", 
						"Start_SELA","Stop_SELA", "Start_CTRB", "Stop_CTRB", "Freq_SELA", "Freq_CTRB", "Pval", "Pval_BF", 
						"Start_SELA","Stop_SELA", "Start_CTRC", "Stop_CTRC", "Freq_SELA", "Freq_CTRC", "Pval", "Pval_BF")
write.table(sa, file=paste(parent.o, "Freqs_SA-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), quote=F, sep="\t", 
							row.names=F, col.names=T)

# Selection pool B
list.files.sel  <- dir(parent.o, pattern ="Freqs_SB-C.*_f0.10.Pval_BF_reads.txt")

input.sa <- read.table(file=paste(parent.o, list.files.sel[1], sep=""), sep="\t", header=T)
input.sb <- read.table(file=paste(parent.o, list.files.sel[2], sep=""), sep="\t", header=T)
input.sc <- read.table(file=paste(parent.o, list.files.sel[3], sep=""), sep="\t", header=T)

sall <- r.sel <- sa <- NULL
sall <- rbind(input.sa, input.sb, input.sc)
sall.u <- sall %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
p.ranges <- GRanges(sall.u$Chr_TE_ID, IRanges(sall.u$Start_SEL, sall.u$Stop_SEL))
p.ranges.red <- reduce(p.ranges)
p.new.int <- as.data.frame(p.ranges.red)
r.sel <- cbind(unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 1)), 
				unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 2)), p.new.int[,2], p.new.int[,3])  
r.sel1 <- as.data.frame(r.sel, stringsAsFactors=FALSE) 
r.sel1[,3] <- as.numeric(r.sel1[,3])
r.sel1[,4] <- as.numeric(r.sel1[,4])
colnames(r.sel1) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new")

sa <- cbind(r.sel1, as.data.frame(matrix("NA", nrow=length(p.new.int[,1]), ncol=24), stringsAsFactors=FALSE))  

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sa[,1])) {
		if (r.sel1[i,1] == input.sa[j,1] && r.sel1[i,2] == input.sa[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sa[j,5], input.sa[j,6])) {
				sa[i, c(5:12)] <- input.sa[j, c(3:10)]
			}
	    }
	}	
}				

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sb[,1])) {
		if (r.sel1[i,1] == input.sb[j,1] && r.sel1[i,2] == input.sb[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sb[j,5], input.sb[j,6])) {
				sa[i, c(13:20)] <- input.sb[j, c(3:10)]
			}
	    }
	}	
}

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sc[,1])) {
		if (r.sel1[i,1] == input.sc[j,1] && r.sel1[i,2] == input.sc[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sc[j,5], input.sc[j,6])) {
				sa[i, c(21:28)] <- input.sc[j, c(3:10)]
			}
	    }
	}	
}

sa[, 12] <- gsub("1", "-", sa[, 12])
sa[, 12] <- gsub("2", "*", sa[, 12])
sa[, 12] <- gsub("3", "**", sa[, 12])
sa[, 20] <- gsub("1", "-", sa[, 20])
sa[, 20] <- gsub("2", "*", sa[, 20])
sa[, 20] <- gsub("3", "**", sa[, 20])
sa[, 28] <- gsub("1", "-", sa[, 28])
sa[, 28] <- gsub("2", "*", sa[, 28])
sa[, 28] <- gsub("3", "**", sa[, 28])

colnames(sa)[5:28] <- c("Start_SELB","Stop_SELB", "Start_CTRA", "Stop_CTRA", "Freq_SELB", "Freq_CTRA", "Pval", "Pval_BF", 
						"Start_SELB","Stop_SELB", "Start_CTRB", "Stop_CTRB", "Freq_SELB", "Freq_CTRB", "Pval", "Pval_BF", 
						"Start_SELB","Stop_SELB", "Start_CTRC", "Stop_CTRC", "Freq_SELB", "Freq_CTRC", "Pval", "Pval_BF")
write.table(sa, file=paste(parent.o, "Freqs_SB-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), quote=F, sep="\t", 
							row.names=F, col.names=T)

# Selection pool C
list.files.sel  <- dir(parent.o, pattern ="Freqs_SC-C.*_f0.10.Pval_BF_reads.txt")

input.sa <- read.table(file=paste(parent.o, list.files.sel[1], sep=""), sep="\t", header=T)
input.sb <- read.table(file=paste(parent.o, list.files.sel[2], sep=""), sep="\t", header=T)
input.sc <- read.table(file=paste(parent.o, list.files.sel[3], sep=""), sep="\t", header=T)

sall <- r.sel <- sa <- NULL
sall <- rbind(input.sa, input.sb, input.sc)
sall.u <- sall %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
p.ranges <- GRanges(sall.u$Chr_TE_ID, IRanges(sall.u$Start_SEL, sall.u$Stop_SEL))
p.ranges.red <- reduce(p.ranges)
p.new.int <- as.data.frame(p.ranges.red)
r.sel <- cbind(unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 1)), 
				unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 2)), p.new.int[,2], p.new.int[,3])  
r.sel1 <- as.data.frame(r.sel, stringsAsFactors=FALSE) 
r.sel1[,3] <- as.numeric(r.sel1[,3])
r.sel1[,4] <- as.numeric(r.sel1[,4])
colnames(r.sel1) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new")

sa <- cbind(r.sel1, as.data.frame(matrix("NA", nrow=length(p.new.int[,1]), ncol=24), stringsAsFactors=FALSE))  

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sa[,1])) {
		if (r.sel1[i,1] == input.sa[j,1] && r.sel1[i,2] == input.sa[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sa[j,5], input.sa[j,6])) {
				sa[i, c(5:12)] <- input.sa[j, c(3:10)]
			}
	    }
	}	
}				

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sb[,1])) {
		if (r.sel1[i,1] == input.sb[j,1] && r.sel1[i,2] == input.sb[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sb[j,5], input.sb[j,6])) {
				sa[i, c(13:20)] <- input.sb[j, c(3:10)]
			}
	    }
	}	
}

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(input.sc[,1])) {
		if (r.sel1[i,1] == input.sc[j,1] && r.sel1[i,2] == input.sc[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(input.sc[j,5], input.sc[j,6])) {
				sa[i, c(21:28)] <- input.sc[j, c(3:10)]
			}
	    }
	}	
}

sa[, 12] <- gsub("1", "-", sa[, 12])
sa[, 12] <- gsub("2", "*", sa[, 12])
sa[, 12] <- gsub("3", "**", sa[, 12])
sa[, 20] <- gsub("1", "-", sa[, 20])
sa[, 20] <- gsub("2", "*", sa[, 20])
sa[, 20] <- gsub("3", "**", sa[, 20])
sa[, 28] <- gsub("1", "-", sa[, 28])
sa[, 28] <- gsub("2", "*", sa[, 28])
sa[, 28] <- gsub("3", "**", sa[, 28])

colnames(sa)[5:28] <- c("Start_SELC","Stop_SELC", "Start_CTRA", "Stop_CTRA", "Freq_SELC", "Freq_CTRA", "Pval", "Pval_BF", 
						"Start_SELC","Stop_SELC", "Start_CTRB", "Stop_CTRB", "Freq_SELC", "Freq_CTRB", "Pval", "Pval_BF", 
						"Start_SELC","Stop_SELC", "Start_CTRC", "Stop_CTRC", "Freq_SELC", "Freq_CTRC", "Pval", "Pval_BF")
write.table(sa, file=paste(parent.o, "Freqs_SC-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), quote=F, sep="\t", 
						row.names=F, col.names=T)

### Classify the TEs of each selection and control pool
# TEs that are significant in three controls and change the frequency in the same direction (increase or decrease)
# TEs that are significant in three controls and change the frequency in the same direction (increase or decrease) and that change the frequency in the same direction in a the third pool but not significantly
# TEs that are significant in one control (increase or decrease) and that change the frequency in the same direction in a the two other pools but not significantly
# All that was different from above

saa <- read.table(file=paste(parent.o, "Freqs_SA-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), header=T, sep="\t")
sbb <- read.table(file=paste(parent.o, "Freqs_SB-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), header=T, sep="\t")
scc <- read.table(file=paste(parent.o, "Freqs_SC-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), header=T, sep="\t")

colnames(saa) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL11","Stop_SEL11","Start_CTR1","Stop_CTR1",
					"Freq_SEL11","Freq_CTR1", "Pval11","Pval_BF11","Start_SEL12","Stop_SEL12","Start_CTR2","Stop_CTR2",
					"Freq_SEL12","Freq_CTR2","Pval12","Pval_BF12","Start_SEL13","Stop_SEL13","Start_CTR3","Stop_CTR3",
					"Freq_SEL13","Freq_CTR3","Pval13","Pval_BF13")
colnames(sbb) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL11","Stop_SEL11","Start_CTR1","Stop_CTR1",
					"Freq_SEL11","Freq_CTR1", "Pval11","Pval_BF11","Start_SEL12","Stop_SEL12","Start_CTR2","Stop_CTR2",
					"Freq_SEL12","Freq_CTR2","Pval12","Pval_BF12","Start_SEL13","Stop_SEL13","Start_CTR3","Stop_CTR3",
					"Freq_SEL13","Freq_CTR3","Pval13","Pval_BF13")
colnames(scc) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL11","Stop_SEL11","Start_CTR1","Stop_CTR1",
					"Freq_SEL11","Freq_CTR1", "Pval11","Pval_BF11","Start_SEL12","Stop_SEL12","Start_CTR2","Stop_CTR2",
					"Freq_SEL12","Freq_CTR2","Pval12","Pval_BF12","Start_SEL13","Stop_SEL13","Start_CTR3","Stop_CTR3",
					"Freq_SEL13","Freq_CTR3","Pval13","Pval_BF13")

# SAA pool
saa.tw.i  <- saa.on.i <- saa.tw.d <- saa.on.d <- saa.all.s <- NULL

# Significant in selection pool and in one, two or three controls and increase in frequency from control to selection pool
saa.tr.i <- saa[grepl("*", saa[,12], fixed = TRUE) & saa[,9] > saa[,10] & grepl("*", saa[,20], fixed = TRUE) & 
				saa[,17] > saa[,18] & grepl("*", saa[,28], fixed = TRUE) & saa[,25] > saa[,26], ]

saa.tw.i.1 <- saa[!grepl("*", saa[,12], fixed = TRUE) & (saa[,9] > saa[,10] | is.na(saa[,9]) | is.na(saa[,10])) & 
				grepl("*", saa[,20], fixed = TRUE) & saa[,17] > saa[,18] & grepl("*", saa[,28], fixed = TRUE) & 
				saa[,25] > saa[,26], ]
saa.tw.i.2 <- saa[grepl("*", saa[,12], fixed = TRUE) & saa[,9] > saa[,10] & !grepl("*", saa[,20], fixed = TRUE) & 
				(saa[,17] > saa[,18] | is.na(saa[,17]) | is.na(saa[,18])) & grepl("*", saa[,28], fixed = TRUE) & 
				saa[,25] > saa[,26], ]
saa.tw.i.3 <- saa[grepl("*", saa[,12], fixed = TRUE) & saa[,9] > saa[,10]  & grepl("*", saa[,20], fixed = TRUE) & 
				saa[,17] > saa[,18] & !grepl("*", saa[,28], fixed = TRUE) & (saa[,25] > saa[,26] | is.na(saa[,25]) | 
				is.na(saa[,26])), ]
saa.tw.i <- rbind(saa.tw.i.1, saa.tw.i.2, saa.tw.i.3)

saa.on.i.1 <- saa[!grepl("*", saa[,12], fixed = TRUE) & (saa[,9] > saa[,10] | is.na(saa[,9]) | is.na(saa[,10])) & 
				!grepl("*", saa[,20], fixed = TRUE) & (saa[,17] > saa[,18] | is.na(saa[,17]) | is.na(saa[,18])) & 
				grepl("*", saa[,28], fixed = TRUE) & saa[,25] > saa[,26], ]
saa.on.i.2 <- saa[!grepl("*", saa[,12], fixed = TRUE) & (saa[,9] > saa[,10] | is.na(saa[,9]) | is.na(saa[,10])) & 
				grepl("*", saa[,20], fixed = TRUE) & saa[,17] > saa[,18] & !grepl("*", saa[,28], fixed = TRUE) & 
				(saa[,25] > saa[,26] | is.na(saa[,25]) | is.na(saa[,26])), ]
saa.on.i.3 <- saa[grepl("*", saa[,12], fixed = TRUE) & saa[,9] > saa[,10] & !grepl("*", saa[,20], fixed = TRUE) & 
				(saa[,17] > saa[,18] | is.na(saa[,17]) | is.na(saa[,18])) & !grepl("*", saa[,28], fixed = TRUE) & 
				(saa[,25] > saa[,26] | is.na(saa[,25]) | is.na(saa[,26])), ]
saa.on.i <- rbind(saa.on.i.1, saa.on.i.2, saa.on.i.3)

# Significant in selection pool and in one, two or three controls and decrease in frequency from control to selection pool
saa.tr.d <- saa[grepl("*", saa[,12], fixed = TRUE) & saa[,9] < saa[,10] & grepl("*", saa[,20], fixed = TRUE) & 
				saa[,17] < saa[,18] & grepl("*", saa[,28], fixed = TRUE) & saa[,25] < saa[,26], ]

saa.tw.d.1 <- saa[!grepl("*", saa[,12], fixed = TRUE) & (saa[,9] < saa[,10] | is.na(saa[,9]) | is.na(saa[,10])) & 
				grepl("*", saa[,20], fixed = TRUE) & saa[,17] < saa[,18] & grepl("*", saa[,28], fixed = TRUE) & 
				saa[,25] < saa[,26], ]
saa.tw.d.2 <- saa[grepl("*", saa[,12], fixed = TRUE) & saa[,9] < saa[,10] & !grepl("*", saa[,20], fixed = TRUE) & 
				(saa[,17] < saa[,18] | is.na(saa[,17]) | is.na(saa[,18])) & grepl("*", saa[,28], fixed = TRUE) & 
				saa[,25] < saa[,26], ]
saa.tw.d.3 <- saa[grepl("*", saa[,12], fixed = TRUE) & saa[,9] < saa[,10] & grepl("*", saa[,20], fixed = TRUE) & 
				saa[,17] < saa[,18] & !grepl("*", saa[,28], fixed = TRUE) & (saa[,25] < saa[,26] | is.na(saa[,25]) | 
				is.na(saa[,26])), ]
saa.tw.d <- rbind(saa.tw.d.1, saa.tw.d.2, saa.tw.d.3)

saa.on.d.1 <- saa[!grepl("*", saa[,12], fixed = TRUE) & (saa[,9] < saa[,10] | is.na(saa[,9]) | is.na(saa[,10])) & 
				!grepl("*", saa[,20], fixed = TRUE) & (saa[,17] < saa[,18] | is.na(saa[,17]) | is.na(saa[,18])) & 
				grepl("*", saa[,28], fixed = TRUE) & saa[,25] < saa[,26], ]
saa.on.d.2 <- saa[!grepl("*", saa[,12], fixed = TRUE) & (saa[,9] < saa[,10] | is.na(saa[,9]) | is.na(saa[,10])) & 
				grepl("*", saa[,20], fixed = TRUE) & saa[,17] < saa[,18] & !grepl("*", saa[,28], fixed = TRUE) & 
				(saa[,25] < saa[,26] | is.na(saa[,25]) | is.na(saa[,26])), ]
saa.on.d.3 <- saa[grepl("*", saa[,12], fixed = TRUE) & saa[,9] < saa[,10] & !grepl("*", saa[,20], fixed = TRUE) & 
				(saa[,17] < saa[,18] | is.na(saa[,17]) | is.na(saa[,18])) & !grepl("*", saa[,28], fixed = TRUE) & 
				(saa[,25] < saa[,26] | is.na(saa[,25]) | is.na(saa[,26])), ]
saa.on.d <- rbind(saa.on.d.1, saa.on.d.2, saa.on.d.3)

# Significant in selection pool and in one, two or three controls and increase/decrease in frequency from control to selection pool
saa.all.s <- rbind(saa.tr.i, saa.tw.i, saa.on.i, saa.tr.d, saa.tw.d, saa.on.d)

# All the rest
saa.o <- anti_join(saa, saa.all.s)

# SBB pool
sbb.tw.i  <- sbb.on.i <- sbb.tw.d <- sbb.on.d <- sbb.all.s <- NULL

# Significant between selection pool and in one, two or three controls and increase in frequency from control to selection pool
sbb.tr.i <- sbb[grepl("*", sbb[,12], fixed = TRUE) & sbb[,9] > sbb[,10] & grepl("*", sbb[,20], fixed = TRUE) & sbb[,17] > sbb[,18] &
				grepl("*", sbb[,28], fixed = TRUE) & sbb[,25] > sbb[,26], ]

sbb.tw.i.1 <- sbb[!grepl("*", sbb[,12], fixed = TRUE) & (sbb[,9] > sbb[,10] | is.na(sbb[,9]) | is.na(sbb[,10])) & 
				grepl("*", sbb[,20], fixed = TRUE) & sbb[,17] > sbb[,18] & grepl("*", sbb[,28], fixed = TRUE) & 
				sbb[,25] > sbb[,26], ]
sbb.tw.i.2 <- sbb[grepl("*", sbb[,12], fixed = TRUE) & sbb[,9] > sbb[,10] & !grepl("*", sbb[,20], fixed = TRUE) & 
				(sbb[,17] > sbb[,18] | is.na(sbb[,17]) | is.na(sbb[,18])) & grepl("*", sbb[,28], fixed = TRUE) & 
				sbb[,25] > sbb[,26], ]
sbb.tw.i.3 <- sbb[grepl("*", sbb[,12], fixed = TRUE) & sbb[,9] > sbb[,10]  & grepl("*", sbb[,20], fixed = TRUE) & 
				sbb[,17] > sbb[,18] & !grepl("*", sbb[,28], fixed = TRUE) & (sbb[,25] > sbb[,26] | is.na(sbb[,25]) | 
				is.na(sbb[,26])), ]
sbb.tw.i <- rbind(sbb.tw.i.1, sbb.tw.i.2, sbb.tw.i.3)

sbb.on.i.1 <- sbb[!grepl("*", sbb[,12], fixed = TRUE) & (sbb[,9] > sbb[,10] | is.na(sbb[,9]) | is.na(sbb[,10])) & 
				!grepl("*", sbb[,20], fixed = TRUE) & (sbb[,17] > sbb[,18] | is.na(sbb[,17]) | is.na(sbb[,18])) & 
				grepl("*", sbb[,28], fixed = TRUE) & sbb[,25] > sbb[,26], ]
sbb.on.i.2 <- sbb[!grepl("*", sbb[,12], fixed = TRUE) & (sbb[,9] > sbb[,10] | is.na(sbb[,9]) | is.na(sbb[,10])) & 
				grepl("*", sbb[,20], fixed = TRUE) & sbb[,17] > sbb[,18] & !grepl("*", sbb[,28], fixed = TRUE) & 
				(sbb[,25] > sbb[,26] | is.na(sbb[,25]) | is.na(sbb[,26])), ]
sbb.on.i.3 <- sbb[grepl("*", sbb[,12], fixed = TRUE) & sbb[,9] > sbb[,10] & !grepl("*", sbb[,20], fixed = TRUE) & 
				(sbb[,17] > sbb[,18] | is.na(sbb[,17]) | is.na(sbb[,18])) & !grepl("*", sbb[,28], fixed = TRUE) & 
				(sbb[,25] > sbb[,26] | is.na(sbb[,25]) | is.na(sbb[,26])), ]
sbb.on.i <- rbind(sbb.on.i.1, sbb.on.i.2, sbb.on.i.3)

# Significant between selection pool and in one, two or three controls and decrease in frequency from control to selection pool
sbb.tr.d <- sbb[grepl("*", sbb[,12], fixed = TRUE) & sbb[,9] < sbb[,10] & grepl("*", sbb[,20], fixed = TRUE) & 
				sbb[,17] < sbb[,18] & grepl("*", sbb[,28], fixed = TRUE) & sbb[,25] < sbb[,26], ]

sbb.tw.d.1 <- sbb[!grepl("*", sbb[,12], fixed = TRUE) & (sbb[,9] < sbb[,10] | is.na(sbb[,9]) | is.na(sbb[,10])) & 
				grepl("*", sbb[,20], fixed = TRUE) & sbb[,17] < sbb[,18] & grepl("*", sbb[,28], fixed = TRUE) & 
				sbb[,25] < sbb[,26], ]
sbb.tw.d.2 <- sbb[grepl("*", sbb[,12], fixed = TRUE) & sbb[,9] < sbb[,10] & !grepl("*", sbb[,20], fixed = TRUE) & 
				(sbb[,17] < sbb[,18] | is.na(sbb[,17]) | is.na(sbb[,18])) & grepl("*", sbb[,28], fixed = TRUE) & 
				sbb[,25] < sbb[,26], ]
sbb.tw.d.3 <- sbb[grepl("*", sbb[,12], fixed = TRUE) & sbb[,9] < sbb[,10] & grepl("*", sbb[,20], fixed = TRUE) & 
				sbb[,17] < sbb[,18] & !grepl("*", sbb[,28], fixed = TRUE) & (sbb[,25] < sbb[,26] | is.na(sbb[,25]) | 
				is.na(sbb[,26])), ]
sbb.tw.d <- rbind(sbb.tw.d.1, sbb.tw.d.2, sbb.tw.d.3)

sbb.on.d.1 <- sbb[!grepl("*", sbb[,12], fixed = TRUE) & (sbb[,9] < sbb[,10] | is.na(sbb[,9]) | is.na(sbb[,10])) & 
				!grepl("*", sbb[,20], fixed = TRUE) & (sbb[,17] < sbb[,18] | is.na(sbb[,17]) | is.na(sbb[,18])) & 
				grepl("*", sbb[,28], fixed = TRUE) & sbb[,25] < sbb[,26], ]
sbb.on.d.2 <- sbb[!grepl("*", sbb[,12], fixed = TRUE) & (sbb[,9] < sbb[,10] | is.na(sbb[,9]) | is.na(sbb[,10])) & 
				grepl("*", sbb[,20], fixed = TRUE) & sbb[,17] < sbb[,18] & !grepl("*", sbb[,28], fixed = TRUE) & 
				(sbb[,25] < sbb[,26] | is.na(sbb[,25]) | is.na(sbb[,26])), ]
sbb.on.d.3 <- sbb[grepl("*", sbb[,12], fixed = TRUE) & sbb[,9] < sbb[,10] & !grepl("*", sbb[,20], fixed = TRUE) & 
				(sbb[,17] < sbb[,18] | is.na(sbb[,17]) | is.na(sbb[,18])) & !grepl("*", sbb[,28], fixed = TRUE) & 
				(sbb[,25] < sbb[,26] | is.na(sbb[,25]) | is.na(sbb[,26])), ]
sbb.on.d <- rbind(sbb.on.d.1, sbb.on.d.2, sbb.on.d.3)

# Significant between selection pool and in one, two or three controls and increase/decrease in frequency from control to selection pool
sbb.all.s <- rbind(sbb.tr.i, sbb.tw.i, sbb.on.i, sbb.tr.d, sbb.tw.d, sbb.on.d)

# All the rest
sbb.o <- anti_join(sbb, sbb.all.s)

# SCC pool
scc.tw.i  <- scc.on.i <- scc.tw.d <- scc.on.d <- scc.all.s <- NULL

# Significant between selection pool and in one, two or three controls and increase in frequency from control to selection pool
scc.tr.i <- scc[grepl("*", scc[,12], fixed = TRUE) & scc[,9] > scc[,10] & grepl("*", scc[,20], fixed = TRUE) & 
				scc[,17] > scc[,18] & grepl("*", scc[,28], fixed = TRUE) & scc[,25] > scc[,26], ]

scc.tw.i.1 <- scc[!grepl("*", scc[,12], fixed = TRUE) & (scc[,9] > scc[,10] | is.na(scc[,9]) | is.na(scc[,10])) & 
				grepl("*", scc[,20], fixed = TRUE) & scc[,17] > scc[,18] & grepl("*", scc[,28], fixed = TRUE) & 
				scc[,25] > scc[,26], ]
scc.tw.i.2 <- scc[grepl("*", scc[,12], fixed = TRUE) & scc[,9] > scc[,10] & !grepl("*", scc[,20], fixed = TRUE) & 
				(scc[,17] > scc[,18] | is.na(scc[,17]) | is.na(scc[,18])) & grepl("*", scc[,28], fixed = TRUE) & 
				scc[,25] > scc[,26], ]
scc.tw.i.3 <- scc[grepl("*", scc[,12], fixed = TRUE) & scc[,9] > scc[,10]  & grepl("*", scc[,20], fixed = TRUE) & 
				scc[,17] > scc[,18] & !grepl("*", scc[,28], fixed = TRUE) & (scc[,25] > scc[,26] | is.na(scc[,25]) | 
				is.na(scc[,26])), ]
scc.tw.i <- rbind(scc.tw.i.1, scc.tw.i.2, scc.tw.i.3)

scc.on.i.1 <- scc[!grepl("*", scc[,12], fixed = TRUE) & (scc[,9] > scc[,10] | is.na(scc[,9]) | is.na(scc[,10])) & 
				!grepl("*", scc[,20], fixed = TRUE) & (scc[,17] > scc[,18] | is.na(scc[,17]) | is.na(scc[,18])) & 
				grepl("*", scc[,28], fixed = TRUE) & scc[,25] > scc[,26], ]
scc.on.i.2 <- scc[!grepl("*", scc[,12], fixed = TRUE) & (scc[,9] > scc[,10] | is.na(scc[,9]) | is.na(scc[,10])) & 
				grepl("*", scc[,20], fixed = TRUE) & scc[,17] > scc[,18] & !grepl("*", scc[,28], fixed = TRUE) & 
				(scc[,25] > scc[,26] | is.na(scc[,25]) | is.na(scc[,26])), ]
scc.on.i.3 <- scc[grepl("*", scc[,12], fixed = TRUE) & scc[,9] > scc[,10] & !grepl("*", scc[,20], fixed = TRUE) & 
				(scc[,17] > scc[,18] | is.na(scc[,17]) | is.na(scc[,18])) & !grepl("*", scc[,28], fixed = TRUE) & 
				(scc[,25] > scc[,26] | is.na(scc[,25]) | is.na(scc[,26])), ]
scc.on.i <- rbind(scc.on.i.1, scc.on.i.2, scc.on.i.3)

# Significant between selection pool and in one, two or three controls and decrease in frequency from control to selection pool
scc.tr.d <- scc[grepl("*", scc[,12], fixed = TRUE) & scc[,9] < scc[,10] & grepl("*", scc[,20], fixed = TRUE) & 
			scc[,17] < scc[,18] & grepl("*", scc[,28], fixed = TRUE) & scc[,25] < scc[,26], ]
			
scc.tw.d.1 <- scc[!grepl("*", scc[,12], fixed = TRUE) & (scc[,9] < scc[,10] | is.na(scc[,9]) | is.na(scc[,10])) & 
				grepl("*", scc[,20], fixed = TRUE) & scc[,17] < scc[,18] & grepl("*", scc[,28], fixed = TRUE) & 
				scc[,25] < scc[,26], ]
scc.tw.d.2 <- scc[grepl("*", scc[,12], fixed = TRUE) & scc[,9] < scc[,10] & !grepl("*", scc[,20], fixed = TRUE) & 
				(scc[,17] < scc[,18] | is.na(scc[,17]) | is.na(scc[,18])) & grepl("*", scc[,28], fixed = TRUE) & 
				scc[,25] < scc[,26], ]
scc.tw.d.3 <- scc[grepl("*", scc[,12], fixed = TRUE) & scc[,9] < scc[,10] & grepl("*", scc[,20], fixed = TRUE) & 
				scc[,17] < scc[,18] & !grepl("*", scc[,28], fixed = TRUE) & (scc[,25] < scc[,26] | is.na(scc[,25]) | 
				is.na(scc[,26])), ]
scc.tw.d <- rbind(scc.tw.d.1, scc.tw.d.2, scc.tw.d.3)

scc.on.d.1 <- scc[!grepl("*", scc[,12], fixed = TRUE) & (scc[,9] < scc[,10] | is.na(scc[,9]) | is.na(scc[,10])) & 
				!grepl("*", scc[,20], fixed = TRUE) & (scc[,17] < scc[,18] | is.na(scc[,17]) | is.na(scc[,18])) & 
				grepl("*", scc[,28], fixed = TRUE) & scc[,25] < scc[,26], ]
scc.on.d.2 <- scc[!grepl("*", scc[,12], fixed = TRUE) & (scc[,9] < scc[,10] | is.na(scc[,9]) | is.na(scc[,10])) & 
				grepl("*", scc[,20], fixed = TRUE) & scc[,17] < scc[,18] & !grepl("*", scc[,28], fixed = TRUE) & 
				(scc[,25] < scc[,26] | is.na(scc[,25]) | is.na(scc[,26])), ]
scc.on.d.3 <- scc[grepl("*", scc[,12], fixed = TRUE) & scc[,9] < scc[,10] & !grepl("*", scc[,20], fixed = TRUE) & 
				(scc[,17] < scc[,18] | is.na(scc[,17]) | is.na(scc[,18])) & !grepl("*", scc[,28], fixed = TRUE) & 
				(scc[,25] < scc[,26] | is.na(scc[,25]) | is.na(scc[,26])), ]
scc.on.d <- rbind(scc.on.d.1, scc.on.d.2, scc.on.d.3)

# Significant between selection pool and in one, two or three controls and increase/decrease in frequency from control to selection pool
scc.all.s <- rbind(scc.tr.i,scc.tw.i, scc.on.i, scc.tr.d, scc.tw.d, scc.on.d)

# All the rest
scc.o <- anti_join(scc, scc.all.s)

# Save all the TEs in three, two or one selections pools ans increase, decrease and all the rest.
# SAA
write.table(rbind(saa.tr.i, saa.tw.i, saa.on.i), file=paste(parent.o, "Freqs_SA-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(rbind(saa.tr.d, saa.tw.d, saa.on.d), file=paste(parent.o, "Freqs_SA-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(saa.o, file=paste(parent.o, "Freqs_SA-3CTRs_Pval_f0.10_BF_reads_others.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)

# SBB
write.table(rbind(sbb.tr.i,sbb.tw.i, sbb.on.i), file=paste(parent.o, "Freqs_SB-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(rbind(sbb.tr.d,sbb.tw.d, sbb.on.d), file=paste(parent.o, "Freqs_SB-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(sbb.o, file=paste(parent.o, "Freqs_SB-3CTRs_f0.10_BF_reads_others.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)

# SCC
write.table(rbind(scc.tr.i,scc.tw.i, scc.on.i), file=paste(parent.o, "Freqs_SC-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(rbind(scc.tr.d,scc.tw.d, scc.on.d), file=paste(parent.o, "Freqs_SC-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(scc.o, file=paste(parent.o, "Freqs_SC-3CTRs_f0.10_BF_reads_others.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)

# Get the TEs that are shared in three, two and one (unique) selection pools (increase and decrease in frequency)
# Increase in frequency

saa <- sbb <- scc <- NULL
saa <- read.table(file=paste(parent.o, "Freqs_SA-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), header=T, sep="\t")
sbb <- read.table(file=paste(parent.o, "Freqs_SB-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), header=T, sep="\t")
scc <- read.table(file=paste(parent.o, "Freqs_SC-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), header=T, sep="\t")

colnames(saa) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR", "Pval","Pval_BF","Start_SEL","Stop_SEL","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR","Pval","Pval_BF","Start_SEL","Stop_SELC","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR","Pval","Pval_BF")
colnames(sbb) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR", "Pval","Pval_BF","Start_SEL","Stop_SEL","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR","Pval","Pval_BF","Start_SEL","Stop_SELC","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR","Pval","Pval_BF")  
colnames(scc) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR", "Pval","Pval_BF","Start_SEL","Stop_SEL","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR","Pval","Pval_BF","Start_SEL","Stop_SELC","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR","Pval","Pval_BF") 

saa[ , 12] <-  gsub("^\\*$", "1", saa[, 12])
saa[ , 12] <-  gsub("^\\**$", "2", saa[, 12])
saa[ , 12] <-  gsub("-", "30", saa[, 12])
saa[ , 12] <- as.numeric(saa[ , 12])
saa[ , 20] <-  gsub("^\\*$", "1", saa[, 20])
saa[ , 20] <-  gsub("^\\**$", "2", saa[, 20])
saa[ , 20] <-  gsub("-", "30", saa[, 20])
saa[ , 20] <- as.numeric(saa[ , 20])
saa[ , 28] <-  gsub("^\\*$", "1", saa[, 28])
saa[ , 28] <-  gsub("^\\**$", "2", saa[, 28])
saa[ , 28] <-  gsub("-", "30", saa[, 28])
saa[ , 28] <- as.numeric(saa[ , 28])

sbb[ , 12] <-  gsub("^\\*$", "1", sbb[, 12])
sbb[ , 12] <-  gsub("^\\**$", "2", sbb[, 12])
sbb[ , 12] <-  gsub("-", "30", sbb[, 12])
sbb[ , 12] <- as.numeric(sbb[ , 12])
sbb[ , 20] <-  gsub("^\\*$", "1", sbb[, 20])
sbb[ , 20] <-  gsub("^\\**$", "2", sbb[, 20])
sbb[ , 20] <-  gsub("-", "30", sbb[, 20])
sbb[ , 20] <- as.numeric(sbb[ , 20])
sbb[ , 28] <-  gsub("^\\*$", "1", sbb[, 28])
sbb[ , 28] <-  gsub("^\\**$", "2", sbb[, 28])
sbb[ , 28] <-  gsub("-", "30", sbb[, 28])
sbb[ , 28] <- as.numeric(sbb[ , 28])
	
scc[ , 12] <-  gsub("^\\*$", "1", scc[, 12])
scc[ , 12] <-  gsub("^\\**$", "2", scc[, 12])
scc[ , 12] <-  gsub("-", "30", scc[, 12])
scc[ , 12] <- as.numeric(scc[ , 12])
scc[ , 20] <-  gsub("^\\*$", "1", scc[, 20])
scc[ , 20] <-  gsub("^\\**$", "2", scc[, 20])
scc[ , 20] <-  gsub("-", "30", scc[, 20])
scc[ , 20] <- as.numeric(scc[ , 20])
scc[ , 28] <-  gsub("^\\*$", "1", scc[, 28])
scc[ , 28] <-  gsub("^\\**$", "2", scc[, 28])
scc[ , 28] <-  gsub("-", "30", scc[, 28])
scc[ , 28] <- as.numeric(scc[ , 28])
 
sall <- NULL
sall <- rbind(saa, sbb, scc)
sall[,1] <- as.character(sall[,1])
sall[,2] <- as.character(sall[,2])
sall <- sall[, 1:4]

sall.u <- sall %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
p.ranges <- GRanges(sall.u$Chr_TE_ID, IRanges(sall.u$Start_SEL_new, sall.u$Stop_SEL_new))
p.ranges.red <- reduce(p.ranges)
p.new.int <- as.data.frame(p.ranges.red)
r.sel <- cbind(unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 1)), 
				unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 2)), p.new.int[,2], p.new.int[,3])  
r.sel1 <- as.data.frame(r.sel, stringsAsFactors=FALSE) 
r.sel1[,3] <- as.numeric(r.sel1[,3])
r.sel1[,4] <- as.numeric(r.sel1[,4])
colnames(r.sel1) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new")

sa <- NULL
sa <- cbind(r.sel1, as.data.frame(matrix("NA", nrow=length(r.sel1[,1]), ncol=78), stringsAsFactors=FALSE))
for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(saa[,1])) {
		if (r.sel1[i,1] == saa[j,1] && r.sel1[i,2] == saa[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(saa[j,3], saa[j,4])) {
				sa[i, c(5:30)] <- (saa[j, c(3:28)])
			}
	    }
	}	
}	

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(sbb[,1])) {
		if (r.sel1[i,1] == sbb[j,1] && r.sel1[i,2] == sbb[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(sbb[j,3], sbb[j,4])) {
				sa[i, c(31:56)] <- (sbb[j, c(3:28)])
			}
	    }
	}	
}

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(scc[,1])) {
		if (r.sel1[i,1] == scc[j,1] && r.sel1[i,2] == scc[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(scc[j,3], scc[j,4])) {
				sa[i, c(57:82)] <- (scc[j, c(3:28)])
			}
	    }
	}	
}			
colnames(sa) <- c("Chr","TE_ID","Start_SEL_all","Stop_SEL_all","Start_SELA_new","Stop_SELA_new","Start_SELA",
				"Stop_SELA","Start_CTRA","Stop_CTRA","Freq_SELA","Freq_CTRA","Pval","Pval_BF","Start_SELA",
				"Stop_SELA","Start_CTRB","Stop_CTRB","Freq_SELA","Freq_CTRB","Pval","Pval_BF","Start_SELA",
				"Stop_SELA","Start_CTRC","Stop_CTRC","Freq_SELA","Freq_CTRC","Pval","Pval_BF","Start_SELB_new",
				"Stop_SELB_new","Start_SELB","Stop_SELB","Start_CTRA","Stop_CTRA","Freq_SELB","Freq_CTRA","Pval",
				"Pval_BF","Start_SELB","Stop_SELB","Start_CTRB","Stop_CTRB","Freq_SELB","Freq_CTRB","Pval","Pval_BF",
				"Start_SELB","Stop_SELB","Start_CTRC","Stop_CTRC","Freq_SELB","Freq_CTRC","Pval","Pval_BF",
				"Start_SELC_new","Stop_SELC_new","Start_SELC","Stop_SELC","Start_CTRA","Stop_CTRA","Freq_SELC",
				"Freq_CTRA","Pval","Pval_BF","Start_SELC","Stop_SELC","Start_CTRB","Stop_CTRB","Freq_SELC",
				"Freq_CTRB","Pval","Pval_BF","Start_SELC","Stop_SELC","Start_CTRC","Stop_CTRC","Freq_SELC",
				"Freq_CTRC","Pval","Pval_BF")

sa[ , 14] <-  gsub("1", "*",  sa[, 14])
sa[ , 14] <-  gsub("2", "**", sa[, 14])
sa[ , 14] <-  gsub("30", "-", sa[, 14])
sa[ , 22] <-  gsub("1", "*", sa[, 22])
sa[ , 22] <-  gsub("2", "**", sa[, 22])
sa[ , 22] <-  gsub("30", "-", sa[, 22])
sa[ , 30] <-  gsub("1", "*", sa[, 30])
sa[ , 30] <-  gsub("2", "**", sa[, 30])
sa[ , 30] <-  gsub("30", "-", sa[, 30])
sa[ , 40] <-  gsub("1", "*", sa[, 40])
sa[ , 40] <-  gsub("2", "**", sa[, 40])
sa[ , 40] <-  gsub("30", "-", sa[, 40])
sa[ , 48] <-  gsub("1", "*", sa[, 48])
sa[ , 48] <-  gsub("2", "**", sa[, 48])
sa[ , 48] <-  gsub("30", "-", sa[, 48])
sa[ , 56] <-  gsub("1", "*", sa[, 56])
sa[ , 56] <-  gsub("2", "**", sa[, 56])
sa[ , 56] <-  gsub("30", "-", sa[, 56])
sa[ , 66] <-  gsub("1", "*", sa[, 66])
sa[ , 66] <-  gsub("2", "**", sa[, 66])
sa[ , 66] <-  gsub("30", "-", sa[, 66])
sa[ , 74] <-  gsub("1", "*", sa[, 74])
sa[ , 74] <-  gsub("2", "**", sa[, 74])
sa[ , 74] <-  gsub("30", "-", sa[, 74])
sa[ , 82] <-  gsub("1", "*", sa[, 82])
sa[ , 82] <-  gsub("2", "**", sa[, 82])
sa[ , 82] <-  gsub("30", "-", sa[, 82])

write.table(sa, file=paste(parent.o, "Freqs_All_SEL-3CTRs_Pval_BF_f0.10_reads_increase.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
sa <- read.table(file=paste(parent.o, "Freqs_All_SEL-3CTRs_Pval_BF_f0.10_reads_increase.txt", sep=""),
				header=T, sep="\t")

on <- tw <- tr <- NULL
for (i in 1:length(sa[,1])) {

	if (!is.na(as.numeric(sa[i,5])) & !is.na(as.numeric(sa[i,31])) & !is.na(as.numeric(sa[i,57]))) {
		tr <- rbind(tr, sa[i, ])
	}
	if (is.na(as.numeric(sa[i,5])) & !is.na(as.numeric(sa[i,31])) & !is.na(as.numeric(sa[i,57]))) {
		tw <- rbind(tw, sa[i, ])
	}
	if (!is.na(as.numeric(sa[i,5])) & !is.na(as.numeric(sa[i,31])) & is.na(as.numeric(sa[i,57]))) {
		tw <- rbind(tw, sa[i, ])
	}
	if (!is.na(as.numeric(sa[i,5])) & is.na(as.numeric(sa[i,31])) & !is.na(as.numeric(sa[i,57]))) {
		tw <- rbind(tw, sa[i, ])
	}
	if (!is.na(as.numeric(sa[i,5])) & is.na(as.numeric(sa[i,31])) & is.na(as.numeric(sa[i,57]))) {
		on <- rbind(on, sa[i, ])
	}
	if (is.na(as.numeric(sa[i,5])) & !is.na(as.numeric(sa[i,31])) & is.na(as.numeric(sa[i,57]))) {
		on <- rbind(on, sa[i, ])
	}
	if (is.na(as.numeric(sa[i,5])) & is.na(as.numeric(sa[i,31])) & !is.na(as.numeric(sa[i,57]))) {
		on <- rbind(on, sa[i, ])
	}
}	
write.table(tr, file=paste(parent.o, "Freqs_3pools_SEL-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(tw, file=paste(parent.o, "Freqs_2pools_SEL-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(on, file=paste(parent.o, "Freqs_1pools_SEL-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)

# Decrease in frequency
saa <- sbb <- scc <- NULL
saa <- read.table(file=paste(parent.o, "Freqs_SA-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), header=T, sep="\t")
sbb <- read.table(file=paste(parent.o, "Freqs_SB-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), header=T, sep="\t")
scc <- read.table(file=paste(parent.o, "Freqs_SC-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), header=T, sep="\t")

colnames(saa) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL","Start_CTR",
					"Stop_CTR","Freq_SEL","Freq_CTR", "Pval","Pval_BF","Start_SEL","Stop_SEL","Start_CTR",   
					"Stop_CTR","Freq_SEL","Freq_CTR","Pval","Pval_BF","Start_SEL","Stop_SELC","Start_CTR",
					"Stop_CTR","Freq_SEL","Freq_CTR","Pval","Pval_BF")
colnames(sbb) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL","Start_CTR",
					"Stop_CTR","Freq_SEL","Freq_CTR", "Pval","Pval_BF","Start_SEL","Stop_SEL","Start_CTR",   
					"Stop_CTR","Freq_SEL","Freq_CTR","Pval","Pval_BF","Start_SEL","Stop_SELC","Start_CTR",
					"Stop_CTR","Freq_SEL","Freq_CTR","Pval","Pval_BF")  
colnames(scc) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR", "Pval","Pval_BF","Start_SEL","Stop_SEL","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR","Pval","Pval_BF","Start_SEL","Stop_SELC","Start_CTR","Stop_CTR",
					"Freq_SEL","Freq_CTR","Pval","Pval_BF") 

saa[ , 12] <-  gsub("^\\*$", "1", saa[, 12])
saa[ , 12] <-  gsub("^\\**$", "2", saa[, 12])
saa[ , 12] <-  gsub("-", "30", saa[, 12])
saa[ , 12] <- as.numeric(saa[ , 12])
saa[ , 20] <-  gsub("^\\*$", "1", saa[, 20])
saa[ , 20] <-  gsub("^\\**$", "2", saa[, 20])
saa[ , 20] <-  gsub("-", "30", saa[, 20])
saa[ , 20] <- as.numeric(saa[ , 20])
saa[ , 28] <-  gsub("^\\*$", "1", saa[, 28])
saa[ , 28] <-  gsub("^\\**$", "2", saa[, 28])
saa[ , 28] <-  gsub("-", "30", saa[, 28])
saa[ , 28] <- as.numeric(saa[ , 28])

sbb[ , 12] <-  gsub("^\\*$", "1", sbb[, 12])
sbb[ , 12] <-  gsub("^\\**$", "2", sbb[, 12])
sbb[ , 12] <-  gsub("-", "30", sbb[, 12])
sbb[ , 12] <- as.numeric(sbb[ , 12])
sbb[ , 20] <-  gsub("^\\*$", "1", sbb[, 20])
sbb[ , 20] <-  gsub("^\\**$", "2", sbb[, 20])
sbb[ , 20] <-  gsub("-", "30", sbb[, 20])
sbb[ , 20] <- as.numeric(sbb[ , 20])
sbb[ , 28] <-  gsub("^\\*$", "1", sbb[, 28])
sbb[ , 28] <-  gsub("^\\**$", "2", sbb[, 28])
sbb[ , 28] <-  gsub("-", "30", sbb[, 28])
sbb[ , 28] <- as.numeric(sbb[ , 28])
	
scc[ , 12] <-  gsub("^\\*$", "1", scc[, 12])
scc[ , 12] <-  gsub("^\\**$", "2", scc[, 12])
scc[ , 12] <-  gsub("-", "30", scc[, 12])
scc[ , 12] <- as.numeric(scc[ , 12])
scc[ , 20] <-  gsub("^\\*$", "1", scc[, 20])
scc[ , 20] <-  gsub("^\\**$", "2", scc[, 20])
scc[ , 20] <-  gsub("-", "30", scc[, 20])
scc[ , 20] <- as.numeric(scc[ , 20])
scc[ , 28] <-  gsub("^\\*$", "1", scc[, 28])
scc[ , 28] <-  gsub("^\\**$", "2", scc[, 28])
scc[ , 28] <-  gsub("-", "30", scc[, 28])
scc[ , 28] <- as.numeric(scc[ , 28])
 
sall <- NULL
sall <- rbind(saa, sbb, scc)
sall[,1] <- as.character(sall[,1])
sall[,2] <- as.character(sall[,2])
sall <- sall[, 1:4]

sall.u <- sall %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&&", remove = TRUE)
p.ranges <- GRanges(sall.u$Chr_TE_ID, IRanges(sall.u$Start_SEL_new, sall.u$Stop_SEL_new))
p.ranges.red <- reduce(p.ranges)
p.new.int <- as.data.frame(p.ranges.red)
r.sel <- cbind(unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 1)), 
				unlist(sapply(strsplit(as.character(p.new.int[,1]),'&&'), "[", 2)), p.new.int[,2], p.new.int[,3])  
r.sel1 <- as.data.frame(r.sel, stringsAsFactors=FALSE) 
r.sel1[,3] <- as.numeric(r.sel1[,3])
r.sel1[,4] <- as.numeric(r.sel1[,4])
colnames(r.sel1) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new")

sa <- NULL
sa <- cbind(r.sel1, as.data.frame(matrix("NA", nrow=length(r.sel1[,1]), ncol=78), stringsAsFactors=FALSE))
for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(saa[,1])) {
		if (r.sel1[i,1] == saa[j,1] && r.sel1[i,2] == saa[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(saa[j,3], saa[j,4])) {
				sa[i, c(5:30)] <- (saa[j, c(3:28)])
			}
	    }
	}	
}	

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(sbb[,1])) {
		if (r.sel1[i,1] == sbb[j,1] && r.sel1[i,2] == sbb[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(sbb[j,3], sbb[j,4])) {
				sa[i, c(31:56)] <- (sbb[j, c(3:28)])
			}
	    }
	}	
}

for (i in 1:length(r.sel1[,1])) {
	for (j in 1:length(scc[,1])) {
		if (r.sel1[i,1] == scc[j,1] && r.sel1[i,2] == scc[j,2]) {
			if (c(r.sel1[i,3],r.sel1[i,4]) %overlaps% c(scc[j,3], scc[j,4])) {
				sa[i, c(57:82)] <- (scc[j, c(3:28)])
			}
	    }
	}	
}			
colnames(sa) <- c("Chr","TE_ID","Start_SEL_all","Stop_SEL_all","Start_SELA_new","Stop_SELA_new","Start_SELA",
					"Stop_SELA","Start_CTRA","Stop_CTRA","Freq_SELA","Freq_CTRA","Pval","Pval_BF","Start_SELA",
					"Stop_SELA","Start_CTRB","Stop_CTRB","Freq_SELA","Freq_CTRB","Pval","Pval_BF","Start_SELA",
					"Stop_SELA","Start_CTRC","Stop_CTRC","Freq_SELA","Freq_CTRC","Pval","Pval_BF","Start_SELB_new",
					"Stop_SELB_new","Start_SELB","Stop_SELB","Start_CTRA","Stop_CTRA","Freq_SELB","Freq_CTRA",
					"Pval","Pval_BF","Start_SELB","Stop_SELB","Start_CTRB","Stop_CTRB","Freq_SELB","Freq_CTRB",
					"Pval","Pval_BF","Start_SELB","Stop_SELB","Start_CTRC","Stop_CTRC","Freq_SELB","Freq_CTRC",
					"Pval","Pval_BF","Start_SELC_new","Stop_SELC_new","Start_SELC","Stop_SELC","Start_CTRA",
					"Stop_CTRA","Freq_SELC","Freq_CTRA","Pval","Pval_BF","Start_SELC","Stop_SELC","Start_CTRB",
					"Stop_CTRB","Freq_SELC","Freq_CTRB","Pval","Pval_BF","Start_SELC","Stop_SELC","Start_CTRC",
					"Stop_CTRC","Freq_SELC","Freq_CTRC","Pval","Pval_BF")

sa[ , 14] <-  gsub("1", "*",  sa[, 14])
sa[ , 14] <-  gsub("2", "**", sa[, 14])
sa[ , 14] <-  gsub("30", "-", sa[, 14])
sa[ , 22] <-  gsub("1", "*", sa[, 22])
sa[ , 22] <-  gsub("2", "**", sa[, 22])
sa[ , 22] <-  gsub("30", "-", sa[, 22])
sa[ , 30] <-  gsub("1", "*", sa[, 30])
sa[ , 30] <-  gsub("2", "**", sa[, 30])
sa[ , 30] <-  gsub("30", "-", sa[, 30])
sa[ , 40] <-  gsub("1", "*", sa[, 40])
sa[ , 40] <-  gsub("2", "**", sa[, 40])
sa[ , 40] <-  gsub("30", "-", sa[, 40])
sa[ , 48] <-  gsub("1", "*", sa[, 48])
sa[ , 48] <-  gsub("2", "**", sa[, 48])
sa[ , 48] <-  gsub("30", "-", sa[, 48])
sa[ , 56] <-  gsub("1", "*", sa[, 56])
sa[ , 56] <-  gsub("2", "**", sa[, 56])
sa[ , 56] <-  gsub("30", "-", sa[, 56])
sa[ , 66] <-  gsub("1", "*", sa[, 66])
sa[ , 66] <-  gsub("2", "**", sa[, 66])
sa[ , 66] <-  gsub("30", "-", sa[, 66])
sa[ , 74] <-  gsub("1", "*", sa[, 74])
sa[ , 74] <-  gsub("2", "**", sa[, 74])
sa[ , 74] <-  gsub("30", "-", sa[, 74])
sa[ , 82] <-  gsub("1", "*", sa[, 82])
sa[ , 82] <-  gsub("2", "**", sa[, 82])
sa[ , 82] <-  gsub("30", "-", sa[, 82])

write.table(sa, file=paste(parent.o, "Freqs_All_SEL-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
sa <- read.table(file=paste(parent.o, "Freqs_All_SEL-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""),
				header=T, sep="\t")

on <- tw <- tr <- NULL
for (i in 1:length(sa[,1])) {

	if (!is.na(as.numeric(sa[i,5])) & !is.na(as.numeric(sa[i,31])) & !is.na(as.numeric(sa[i,57]))) {
		tr <- rbind(tr, sa[i, ])
	}
	if (is.na(as.numeric(sa[i,5])) & !is.na(as.numeric(sa[i,31])) & !is.na(as.numeric(sa[i,57]))) {
		tw <- rbind(tw, sa[i, ])
	}
	if (!is.na(as.numeric(sa[i,5])) & !is.na(as.numeric(sa[i,31])) & is.na(as.numeric(sa[i,57]))) {
		tw <- rbind(tw, sa[i, ])
	}
	if (!is.na(as.numeric(sa[i,5])) & is.na(as.numeric(sa[i,31])) & !is.na(as.numeric(sa[i,57]))) {
		tw <- rbind(tw, sa[i, ])
	}
	if (!is.na(as.numeric(sa[i,5])) & is.na(as.numeric(sa[i,31])) & is.na(as.numeric(sa[i,57]))) {
		on <- rbind(on, sa[i, ])
	}
	if (is.na(as.numeric(sa[i,5])) & !is.na(as.numeric(sa[i,31])) & is.na(as.numeric(sa[i,57]))) {
		on <- rbind(on, sa[i, ])
	}
	if (is.na(as.numeric(sa[i,5])) & is.na(as.numeric(sa[i,31])) & !is.na(as.numeric(sa[i,57]))) {
		on <- rbind(on, sa[i, ])
	}
}	
write.table(tr, file=paste(parent.o, "Freqs_3pools_SEL-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(tw, file=paste(parent.o, "Freqs_2pools_SEL-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)
write.table(on, file=paste(parent.o, "Freqs_1pools_SEL-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T)

# Get the frequencies of other selection pools corresponfing to TEs that are significant in one selection pool
saa <- read.table(file=paste(parent.o, "Freqs_SA-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE)
sbb <- read.table(file=paste(parent.o, "Freqs_SB-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE)
scc <- read.table(file=paste(parent.o, "Freqs_SC-3CTRs_f0.10.Pval_BF_reads.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE)

oni <- read.table(file=paste(parent.o, "Freqs_1pools_SEL-3CTRs_Pval_f0.10_BF_reads_increase.txt", sep=""), sep="\t", 
					header=T, stringsAsFactors = FALSE)
ond <- read.table(file=paste(parent.o, "Freqs_1pools_SEL-3CTRs_Pval_f0.10_BF_reads_decrease.txt", sep=""), sep="\t", 
					header=T, stringsAsFactors = FALSE)

sa.i <- oni[!is.na(oni[,6]) & is.na(oni[,31]) & is.na(oni[,57]), ]
sb.i <- oni[is.na(oni[,6]) & !is.na(oni[,31]) & is.na(oni[,57]), ]
sc.i <- oni[is.na(oni[,6]) & is.na(oni[,31]) & !is.na(oni[,57]), ]

sa.d <- ond[!is.na(ond[,6]) & is.na(ond[,31]) & is.na(ond[,57]), ]
sb.d <- ond[is.na(ond[,6]) & !is.na(ond[,31]) & is.na(ond[,57]), ]
sc.d <- ond[is.na(ond[,6]) & is.na(ond[,31]) & !is.na(ond[,57]), ]

## 1D pool
# Increase in frequency
same <- pool <- NULL
for (i in 1:length(sbb[,1])) {
	for (j in 1:length(sa.i[,1])) {
		if (sbb[i,1] == sa.i[j,1] && sbb[i,2] == sa.i[j,2]) {
			if (c(sbb[i,3],sbb[i,4]) %overlaps% c(sa.i[j,3], sa.i[j,4])) {
				same <- rbind(same, sbb[i,])
				pool <- c(pool, "2D")
			} 
		}	
	}	
}

colnames(same) <- colnames(scc)

for (i in 1:length(scc[,1])) {
	for (j in 1:length(sa.i[,1])) {
		if (scc[i,1] == sa.i[j,1] && scc[i,2] == sa.i[j,2]) {
			if (c(scc[i,3],scc[i,4]) %overlaps% c(sa.i[j,3], sa.i[j,4])) {
				same <- rbind(same, scc[i,])
				pool <- c(pool, "3D")
			} 
		}	
	}	
}
saa.i <- cbind(same, pool)
colnames(saa.i) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL",
					"Start_CTRA","Stop_CTRA","Freq_SEL","Freq_CTRA","Pval","Pval_BF",
					"Start_SEL.1","Stop_SEL.1","Start_CTRB","Stop_CTRB","Freq_SEL.1","Freq_CTR", 
					"Pval.1","Pval_BF.1","Start_SEL.2","Stop_SEL.2","Start_CTRC","Stop_CTRC", 
					"Freq_SEL.2","Freq_CTRC","Pval.2","Pval_BF.2", "other_pool_non_sign")

write.table(saa.i, file=paste(parent.o, "One_sign_SA_increase_freqs_in_other_pools_f0.10.Pval_BF_reads.txt", sep=""), 
			quote=F, col.names=T, row.names=F, sep="\t")

# Decrease in frequency
same.d <- pool <- NULL

for (i in 1:length(sbb[,1])) {
	for (j in 1:length(sa.d[,1])) {
		if (sbb[i,1] == sa.d[j,1] && sbb[i,2] == sa.d[j,2]) {
			if (c(sbb[i,3],sbb[i,4]) %overlaps% c(sa.d[j,3], sa.d[j,4])) {
				same.d <- rbind(same.d, sbb[i,])
				pool <- c(pool, "2D")
			} 
		}	
	}	
}

colnames(same.d) <- colnames(scc)

for (i in 1:length(scc[,1])) {
	for (j in 1:length(sa.d[,1])) {
		if (scc[i,1] == sa.d[j,1] && scc[i,2] == sa.d[j,2]) {
			if (c(scc[i,3],scc[i,4]) %overlaps% c(sa.d[j,3], sa.d[j,4])) {
				same.d <- rbind(same.d, scc[i,])
				pool <- c(pool, "3D")
			} 
		}	
	}	
}
saa.d <- cbind(same.d, pool)
colnames(saa.d) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL",
					"Start_CTRA","Stop_CTRA","Freq_SEL","Freq_CTRA","Pval","Pval_BF",
					"Start_SEL.1","Stop_SEL.1","Start_CTRB","Stop_CTRB","Freq_SEL.1","Freq_CTR", 
					"Pval.1","Pval_BF.1","Start_SEL.2","Stop_SEL.2","Start_CTRC","Stop_CTRC", 
					"Freq_SEL.2","Freq_CTRC","Pval.2","Pval_BF.2", "other_pool_non_sign")
write.table(saa.d, file=paste(parent.o, "One_sign_SA_decrease_freqs_in_other_pools_f0.10.Pval_BF_reads.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")

## 2D pool
# Increase in frequency
same <- pool <- NULL
for (i in 1:length(saa[,1])) {
	for (j in 1:length(sb.i[,1])) {
		if (saa[i,1] == sb.i[j,1] && saa[i,2] == sb.i[j,2]) {
			if (c(saa[i,3],saa[i,4]) %overlaps% c(sb.i[j,3], sb.i[j,4])) {
				same <- rbind(same, saa[i,])
				pool <- c(pool, "1D")
			} 
		}	
	}	
}

colnames(same) <- colnames(scc)

for (i in 1:length(scc[,1])) {
	for (j in 1:length(sb.i[,1])) {
		if (scc[i,1] == sb.i[j,1] && scc[i,2] == sb.i[j,2]) {
			if (c(scc[i,3], scc[i,4]) %overlaps% c(sb.i[j,3], sb.i[j,4])) {
				same <- rbind(same, scc[i,])
				pool <- c(pool, "3D")
			} 
		}	
	}	
}

sbb.i <- cbind(same, pool)
colnames(sbb.i) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL",
					"Start_CTRA","Stop_CTRA","Freq_SEL","Freq_CTRA","Pval","Pval_BF",
					"Start_SEL.1","Stop_SEL.1","Start_CTRB","Stop_CTRB","Freq_SEL.1","Freq_CTR", 
					"Pval.1","Pval_BF.1","Start_SEL.2","Stop_SEL.2","Start_CTRC","Stop_CTRC", 
					"Freq_SEL.2","Freq_CTRC","Pval.2","Pval_BF.2", "other_pool_non_sign")

write.table(sbb.i, file=paste(parent.o, "One_sign_SB_increase_freqs_in_other_pools_f0.10.Pval_BF_reads.txt", sep=""), 
			quote=F, col.names=T, row.names=F, sep="\t")

# Decrease in frequency
same.d <- pool <- NULL

for (i in 1:length(saa[,1])) {
	for (j in 1:length(sb.d[,1])) {
		if (saa[i,1] == sb.d[j,1] && saa[i,2] == sb.d[j,2]) {
			if (c(saa[i,3],saa[i,4]) %overlaps% c(sb.d[j,3], sb.d[j,4])) {
				same.d <- rbind(same.d, saa[i,])
				pool <- c(pool, "1D")
			} 
		}	
	}	
}

colnames(same.d) <- colnames(scc)

for (i in 1:length(scc[,1])) {
	for (j in 1:length(sb.d[,1])) {
		if (scc[i,1] == sb.d[j,1] && scc[i,2] == sb.d[j,2]) {
			if (c(scc[i,3],scc[i,4]) %overlaps% c(sb.d[j,3], sb.d[j,4])) {
				same.d <- rbind(same.d, scc[i,])
				pool <- c(pool, "3D")
			} 
		}	
	}	
}
sbb.d <- cbind(same.d, pool)
colnames(sbb.d) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL",
					"Start_CTRA","Stop_CTRA","Freq_SEL","Freq_CTRA","Pval","Pval_BF",
					"Start_SEL.1","Stop_SEL.1","Start_CTRB","Stop_CTRB","Freq_SEL.1","Freq_CTR", 
					"Pval.1","Pval_BF.1","Start_SEL.2","Stop_SEL.2","Start_CTRC","Stop_CTRC", 
					"Freq_SEL.2","Freq_CTRC","Pval.2","Pval_BF.2", "other_pool_non_sign")
write.table(sbb.d, file=paste(parent.o, "One_sign_SB_decrease_freqs_in_other_pools_f0.10.Pval_BF_reads.txt", sep=""), 
			quote=F, col.names=T, row.names=F, sep="\t")

## 3D pool
# Increase in frequency
same <- pool <- NULL
for (i in 1:length(saa[,1])) {
	for (j in 1:length(sc.i[,1])) {
		if (saa[i,1] == sc.i[j,1] && saa[i,2] == sc.i[j,2]) {
			if (c(saa[i,3],saa[i,4]) %overlaps% c(sc.i[j,3], sc.i[j,4])) {
				same <- rbind(same, saa[i,])
				pool <- c(pool, "1D")
			} 
		}	
	}	
}

colnames(same) <- colnames(sbb)

for (i in 1:length(sbb[,1])) {
	for (j in 1:length(sc.i[,1])) {
		if (sbb[i,1] == sc.i[j,1] && sbb[i,2] == sc.i[j,2]) {
			if (c(sbb[i,3],sbb[i,4]) %overlaps% c(sc.i[j,3], sc.i[j,4])) {
				same <- rbind(same, sbb[i,])
				pool <- c(pool, "2D")
			} 
		}	
	}	
}
scc.i <- cbind(same, pool)
colnames(scc.i) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL",
					"Start_CTRA","Stop_CTRA","Freq_SEL","Freq_CTRA","Pval","Pval_BF",
					"Start_SEL.1","Stop_SEL.1","Start_CTRB","Stop_CTRB","Freq_SEL.1","Freq_CTR", 
					"Pval.1","Pval_BF.1","Start_SEL.2","Stop_SEL.2","Start_CTRC","Stop_CTRC", 
					"Freq_SEL.2","Freq_CTRC","Pval.2","Pval_BF.2", "other_pool_non_sign")

write.table(scc.i, file=paste(parent.o, "One_sign_S3_increase_freqs_in_other_pools_f0.10.Pval_BF_reads.txt", sep=""), 
			quote=F, col.names=T, row.names=F, sep="\t")

# Decrease in frequency
same.d <- pool <- NULL

for (i in 1:length(saa[,1])) {
	for (j in 1:length(sc.d[,1])) {
		if (saa[i,1] == sc.d[j,1] && saa[i,2] == sc.d[j,2]) {
			if (c(saa[i,3],saa[i,4]) %overlaps% c(sc.d[j,3], sc.d[j,4])) {
				same.d <- rbind(same.d, saa[i,])
				pool <- c(pool, "1D")
			} 
		}	
	}	
}

colnames(same.d) <- colnames(sbb)

for (i in 1:length(sbb[,1])) {
	for (j in 1:length(sc.d[,1])) {
		if (sbb[i,1] == sc.d[j,1] && sbb[i,2] == sc.d[j,2]) {
			if (c(sbb[i,3],sbb[i,4]) %overlaps% c(sc.d[j,3], sc.d[j,4])) {
				same.d <- rbind(same.d, sbb[i,])
				pool <- c(pool, "2D")
			} 
		}	
	}	
}
scc.d <- cbind(same.d, pool)
colnames(scc.d) <- c("Chr","TE_ID","Start_SEL_new","Stop_SEL_new","Start_SEL","Stop_SEL",
					"Start_CTRA","Stop_CTRA","Freq_SEL","Freq_CTRA","Pval","Pval_BF",
					"Start_SEL.1","Stop_SEL.1","Start_CTRB","Stop_CTRB","Freq_SEL.1","Freq_CTR", 
					"Pval.1","Pval_BF.1","Start_SEL.2","Stop_SEL.2","Start_CTRC","Stop_CTRC", 
					"Freq_SEL.2","Freq_CTRC","Pval.2","Pval_BF.2", "other_pool_non_sign")
write.table(scc.d, file=paste(parent.o, "One_sign_S3_decrease_freqs_in_other_pools_f0.10.Pval_BF_reads.txt", sep=""), quote=F, col.names=T, row.names=F, sep="\t")
