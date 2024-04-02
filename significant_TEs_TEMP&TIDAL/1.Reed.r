#!/usr/bin/env Rscript

library(DescTools) 
library(GenomicRanges)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with the TEMP output files (i.e., "/home/")
#args[2]: path to the results folder 
#args[3]: file with the TEs dictionary FlyBase-Repbase (with full path) 
#args[4]: path to the folder containing the TIDAL output files (i.e., "/home/")

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
	ins.1p1.0.10.ref.1 <- NULL
	for (i in 1:length(ins.1p1.0.10[,1])) {
		if (as.numeric(ins.1p1.0.10[i,10]) != 0 && as.numeric(ins.1p1.0.10[i,12]) != 0 && 
			as.numeric(ins.1p1.0.10[i,13]) != 0 && as.numeric(ins.1p1.0.10[i,14]) != 0) {		
			ins.1p1.0.10.ref.1 <- rbind(ins.1p1.0.10.ref.1, ins.1p1.0.10[i, c(1,9,11,4:8,10,12:14)])		
		} 
	}
		
   # Predicted junctions with bp-resolution (columns 10,11:13 with some value = 0) = Broad junctions
	ins.1p1.0.10.ref.2 <- NULL
	for (i in 1:length(ins.1p1.0.10[,1])) {	
		if (as.numeric(ins.1p1.0.10[i,10]) == 0 || as.numeric(ins.1p1.0.10[i,12]) == 0 || 
			as.numeric(ins.1p1.0.10[i,13]) == 0 || as.numeric(ins.1p1.0.10[i,14]) == 0) {	 
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
	data.f0.10.final <- cbind(matrix(unlist(strsplit(data.f0.10.final.j$All, split="//", fixed=T)), ncol=4, byrow=T, 
								dimnames=list(NULL, c("Chr", "Start", "End", "TransposonName"))), data.f0.10.final.j[, c(2:10)])
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
				y[l,4] <- as.character(unique(nom[as.character(nom[ ,2]) == as.character(y[l,4]), 1])[1]) 
			} else { y[l,4] <- y[l,4] }
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
		write.table(m.final, paste(parent.O, nam, ".final.txt", sep=""), 
					quote=F, col.names=T, row.names=F, sep="\t")
}	

### Read TIDAL output files
temp.folder <- "TEMP/insertion/"
tidal.folder <- args[4]

list.files.temp  <- dir(parent.o, pattern ="*insertion.1p1.f0.10.final.txt")
list.files.tidal  <- dir(tidal.folder, pattern ="*Inserts_Annotated.txt")

for (f in 1:length(list.files.temp)) {
	input.temp <- paste(parent.o, list.files.temp[f], sep="")
	input.tidal <- paste(parent.o, list.files.tidal[f], sep="")
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
						data1.p <- rbind(data1.p, c(as.character(temp[i,1]), min.range, max.range, 
										as.character(temp[i,4]), temp[i,7], temp[i,8]))
					}
				}		
			}
		}
		colnames(data1.p) <- c("Chr", "Start", "Stop", "TE_ID", "Reads", "Freq")
	
		# Remove redundancy
		data1.j <- as.data.frame(data1.p) %>% unite(All, c(Chr, Start,  Stop, TE_ID, Reads, Freq), sep = "&", remove = TRUE)
		data1.u <- unique(data1.j)
		data1 <- matrix(unlist(strsplit(data1.u$All, split="&", fixed=T)), ncol=6, byrow=T, 
						dimnames=list(NULL, c("Chr", "Start", "Stop", "TE_ID", "Reads", "Freq")))

		# Get the number of reads for the reference alleles
		data1n <- data1
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
				mr1 <- mr[as.numeric(as.character(mr[,6])) + as.numeric(as.character(mr[,7])) == max(as.numeric(as.character(mr[,6])) + as.numeric(as.character(mr[,7]))) , ]
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

### Compare BL and SELECTION pools to get the frequencies of the intersected intervals and perform Fisher's exact test and Bonferroni correction
# Build single intervals file for the BL from the 7 different BL pools
list.files.bl.1  <- dir(parent.o, pattern ="BL.*_new_intervals_1p1.f0.10_reads.txt")

# Append the intervals for all BL pools
all.bs.1 <-  NULL
for (i in 1:length(list.files.bl.1)) {
	bs.1 <- read.table(file=paste(parent.o,list.files.bl.1[i], sep=""), header=T, sep="\t") 
	all.bs.1 <- rbind(all.bs.1, bs.1)
}

# Get new intervals from the appended BL pools to remove overlapping.
all.bs.1u <- all.bs.1 %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&", remove = TRUE)
all.bs.1r <- GRanges(all.bs.1u$Chr_TE_ID, IRanges(all.bs.1u$Start, all.bs.1u$Stop))
all.bs.1r.red <- reduce(all.bs.1r)
all.1 <- as.data.frame(all.bs.1r.red)
all.1.s <- all.1[,1:3]
colnames(all.1.s) <- c("Chr_TE_ID", "Start", "Stop")

# Compare the new intervals that do not overlap with each BL pool and get the interval among the BL pools with maximum reads 
dat <- r.bl <- NULL
for (i in 1:length(all.1.s[,1])) {
	dat <- NULL
	for (j in 1:length(list.files.bl.1)) {
		x <- read.table(file=paste(parent.o,list.files.bl.1[j], sep=""), header=T, sep="\t")
		x.u <- x %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "&", remove = TRUE)
		for (k in 1:length(x.u[,1])) { 
			if (all.1.s[i,1] == x.u[k,1]) {
				if (c(all.1.s[i,2], all.1.s[i,3]) %overlaps% c(x.u[k,2], x.u[k,3])) {
					dat <- rbind(dat, x.u[k, ])
				}	
			}				
		}
	}
	if (!is.null(dat) && dim(dat)[1] > 1) { 
		m <- max(rowSums(dat[,5:6]))
		su <- rowSums(dat[,5:6])
		dat.max <- dat[names(su[su == m]), ] 
		if ( dim(dat.max)[1] > 1 ) { 
		 	if (dim(dat.max[dat.max[,4] == max(dat.max[,4]), ])[1] == 1) { 
		 		r.bl <- rbind(r.bl, dat.max[dat.max[,4] == max(dat.max[,4]), ])
		 	} else {
		 		r.bl <- rbind(r.bl,  dat.max[dat.max[,4] == max(dat.max[,4]), ][1,])
		 	}
		} else {
			r.bl <- rbind(r.bl, dat.max)
		}
	}
	if ( !is.null(dat) && dim(dat)[1] == 1) {
		r.bl <- rbind(r.bl, dat)
	}	 	
}	

r.bl.final <- cbind(matrix(unlist(strsplit(as.character(r.bl$Chr_TE_ID), split="&", fixed=T)), 
					ncol=2, byrow=T, dimnames=list(NULL, c("Chr", "TE_ID"))),r.bl[,2:6]) 
write.table(all.1.s, file=paste(parent.o, "concatenated.all.bl.1p1.f0.10.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T) 
write.table(r.bl.final, file=paste(parent.o, "concatenated.all.bl.1p1.f0.10.final.txt", sep=""), 
			quote=F, sep="\t", row.names=F, col.names=T) 

# Compare BL with each SELECTION pool to get the frequencies and perform Fisher's exact test and Bonferroni correction
input.c <- read.table(file=paste(parent.o, "concatenated.all.bl.1p1.f0.10.final.txt", sep=""), sep="\t", header=T)

# HS pools
for (f in 1:2) {	
	input.r <- read.table(file=paste(parent.o, "HS",f,"_new_intervals_1p1.f0.10_reads.txt", sep=""), sep="\t", header=T)
	data1 <- NULL
	for (i in 1:length(input.c[,1])) {
		for (j in 1:length(input.r[,1])) {
			if (input.c[i,1] == input.r[j,1] && input.c[i,2] == input.r[j,4]) {
				if (c(input.c[i,3],input.c[i,4]) %overlaps% c(input.r[j,2], input.r[j,3])) {
					pval <- fisher.test(matrix(c(as.numeric(input.c[i,6]), as.numeric(input.c[i,7]), 
										as.numeric(input.r[j,6]), as.numeric(input.r[j,7])), nrow=2, ncol=2, byrow=T))$p.val
					data1 <- rbind(data1, c(as.character(input.c[i,1]), as.character(input.c[i,2]), as.numeric(input.c[i,3]), 
									as.numeric(input.c[i,4]), as.numeric(input.r[j,2]), as.numeric(input.r[j,3]), as.numeric(input.c[i,5]), 
									as.numeric(input.r[j,5]), scientific(as.numeric(pval), digits=2)))
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
		colnames(data2) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", 
							"Pval", "Pval_BF")
		write.table(data2, file=paste(parent.o, "Freqs_BL-HS", f,"_f0.10_Pval_BF_reads.txt", sep=""), 
										quote=F, sep="\t", row.names=F, col.names=T)
	}
}

# LS pool
for (f in 1:2) {	
	input.r <- read.table(file=paste(parent.o, "LS",f,"_new_intervals_1p1.f0.10_reads.txt", sep=""), sep="\t", header=T)
	data1 <- NULL
	for (i in 1:length(input.c[,1])) {
		for (j in 1:length(input.r[,1])) {
			if (input.c[i,1] == input.r[j,1] && input.c[i,2] == input.r[j,4]) {
				if (c(input.c[i,3],input.c[i,4]) %overlaps% c(input.r[j,2], input.r[j,3])) {
					pval <- fisher.test(matrix(c(as.numeric(input.c[i,6]), as.numeric(input.c[i,7]), 
										as.numeric(input.r[j,6]), as.numeric(input.r[j,7])), nrow=2, ncol=2, byrow=T))$p.val
					data1 <- rbind(data1, c(as.character(input.c[i,1]), as.character(input.c[i,2]), 
									as.numeric(input.c[i,3]), as.numeric(input.c[i,4]), as.numeric(input.r[j,2]), 
									as.numeric(input.r[j,3]), as.numeric(input.c[i,5]), as.numeric(input.r[j,5]), 
									scientific(as.numeric(pval), digits=2)))
				}
			}		
		}
	}
	colnames(data1) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval")
	write.table(data1, file=paste(parent.o, "Freqs_BL-LS", f,"_f0.10_Pval_reads.txt", sep=""), 
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
		colnames(data2) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", 
							"Pval", "Pval_BF")
		write.table(data2, file=paste(parent.o, "Freqs_BL-LS", f,"_f0.10_Pval_BF_reads.txt", sep=""), 
					quote=F, sep="\t", row.names=F, col.names=T)
	}
}

# HF
for (f in 1:2) {	
	input.r <- read.table(file=paste(parent.o, "HF",f,"_new_intervals_1p1.f0.10_reads.txt", sep=""), sep="\t", header=T)
	data1 <- NULL
	for (i in 1:length(input.c[,1])) {
		for (j in 1:length(input.r[,1])) {
			if (input.c[i,1] == input.r[j,1] && input.c[i,2] == input.r[j,4]) {
				if (c(input.c[i,3],input.c[i,4]) %overlaps% c(input.r[j,2], input.r[j,3])) {
					pval <- fisher.test(matrix(c(as.numeric(input.c[i,6]), as.numeric(input.c[i,7]), as.numeric(input.r[j,6]), 
										as.numeric(input.r[j,7])), nrow=2, ncol=2, byrow=T))$p.val
					data1 <- rbind(data1, c(as.character(input.c[i,1]), as.character(input.c[i,2]), as.numeric(input.c[i,3]), 
									as.numeric(input.c[i,4]), as.numeric(input.r[j,2]), as.numeric(input.r[j,3]), 
									as.numeric(input.c[i,5]), as.numeric(input.r[j,5]), scientific(as.numeric(pval), digits=2)))
				}
			}		
		}
	}
	colnames(data1) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval")
	write.table(data1, file=paste(parent.o, "Freqs_BL-HF", f,"_f0.10_Pval_reads.txt", sep=""), quote=F, sep="\t", 
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
		colnames(data2) <- c("Chr", "TE_ID", "Start_BL", "Stop_BL", "Start_S", "Stop_S", "Freq_BL", "Freq_S", "Pval", "Pval_BF")
		write.table(data2, file=paste(parent.o, "Freqs_BL-HF", f,"_f0.10_Pval_BF_reads.txt", sep=""), quote=F, 
					sep="\t", row.names=F, col.names=T)
	}
}
