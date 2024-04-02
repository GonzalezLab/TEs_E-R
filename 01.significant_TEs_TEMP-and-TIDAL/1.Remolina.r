#!/usr/bin/env Rscript

library(DescTools) 
library(GenomicRanges)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with the diferent subfolders (pools) containing the TEMP output files (i.e., "/home/")
#args[2]: file with the TEs dictionary FlyBase-Repbase (with full path)
#args[3]: path to the folder with the diferent subfolders (pools) containing the TIDAL output files (i.e., "/home/")
#args[4]: path to the results folder

### Read TEMP output files
parent.f <- args[1]

list.dirs  <- dir(parent.f, pattern ="*-deSim-mel")

for (d in 1:length(list.dirs)) {
	parent.o <- paste(parent.f, list.dirs[d], "/", sep="")
	setwd(paste(parent.f, list.dirs[d], sep=""))
	input <- list.files(pattern = "summary")
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
		
	# Include the name of the strain for each input file (ex: CSD1, CSD2...RSD1)	
	strain <- rep(name, length(ins.1p1.0.10.ref.1[,1]))
	ins.1p1.0.10.ref.1 <- cbind(ins.1p1.0.10.ref.1, strain)
	strain <- rep(name, length(ins.1p1.0.10.ref.2[,1]))
	ins.1p1.0.10.ref.2 <- cbind(ins.1p1.0.10.ref.2, strain)	

	# Join both the Narrow and Broad junction and name the final data.frame
	ins.1p1.0.10.ref.1 <- unname(ins.1p1.0.10.ref.1)
	ins.1p1.0.10.ref.2 <- unname(ins.1p1.0.10.ref.2)
	colnames(ins.1p1.0.10.ref.1) <- c("Chr", "Start", "End", "TransposonName", "TransposonDirection", "Class", "VariantSupport", "Frequency", "Junction1Support", "Junction2Support", "5'_Support", "3'_Support", "Strain")
	colnames(ins.1p1.0.10.ref.2) <- c("Chr", "Start", "End", "TransposonName", "TransposonDirection", "Class", "VariantSupport", "Frequency", "Junction1Support", "Junction2Support", "5'_Support", "3'_Support", "Strain")
	ins.1p1.0.10.ref <- rbind(ins.1p1.0.10.ref.1, ins.1p1.0.10.ref.2)
	
	# Remove those TEs that are the same in the same chr and coordinates but annotated twice
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
	data.f0.10.final <- cbind(matrix(unlist(strsplit(data.f0.10.final.j$All, split="//", fixed=T)), ncol=4, byrow=T, 
					dimnames=list(NULL, c("Chr", "Start", "End", "TransposonName"))), data.f0.10.final.j[, c(2:10)])
	data.f0.10.final[, 2] <- as.numeric(as.character(data.f0.10.final[,2]))
	data.f0.10.final[, 3] <- as.numeric(as.character(data.f0.10.final[,3]))
		 
	write.table(data.f0.10.final, file=paste(parent.o, output, sep=""), quote=F, col.names=T, row.names=F, sep="\t")	 	
}

### Change nomenclature in TEMP (REPBASE nomenclature for TEs) files for uniformity with TIDAL (FlyBase nomenclature)
parent.f <- args[1]
list.dirs  <- dir(parent.f, pattern ="*-deSim-mel")

nom.file <- args[2]
nom <- read.table(nom.file, header=T, sep="\t")

for (d in 1:length(list.dirs)) {
 	list.dirs.f <- dir(paste(parent.f,list.dirs[d], sep=""), pattern ="*insertion.1p1.*txt")
 	for (f in 1:length(list.dirs.f)) {
 		x <- read.table(file=paste(parent.f,list.dirs[d], "/", list.dirs.f[f], sep=""), 
 						header=T, sep="\t", fill=T, row.names=NULL)
 		colnames(x) <- c("Chr", "Start", "End", "TransposonName", "TransposonDirection", "Class", 
 						"VariantSupport", "Frequency", "Junction1Support", "Junction2Support", 
 						"p5_Support", "p3_Support", "Strain")
 		y <- x
		y[,4] <- as.character(y[,4])
		for (l in 1:length(y[,1])) {
			if (as.character(y[l,4]) %in% as.character(nom[,2])) {
				y[l,4] <- as.character(unique(nom[as.character(nom[ ,2]) == as.character(y[l,4]), 1])[1]) 
			} else { y[l,4] <- y[l,4] }
 		}
 		nam <-strsplit(list.dirs.f[f], split=".txt", fixed=T)[[1]][1]
 
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
		write.table(m.final, paste(parent.f, list.dirs[d], "/", nam, ".final.txt", sep=""), 
					quote=F, col.names=T, row.names=F, sep="\t")
	}	
}

### Read TIDAL output files
tidal.folder <- args[3]
parent.r <- args[4]

list.dirs.temp  <- dir(parent.f, pattern ="*-deSim-mel")
list.dirs.tidal <- dir(tidal.folder, pattern ="*")

# Get new intervals from those for which temp and tidal intersect: 
for (f in 1:length(list.dirs.temp)) {
	input.temp <- list.files(paste(parent.f, list.dirs.temp[f], sep=""), 
					pattern = "*insertion.1p1.f0.10.final.txt")
	input.tidal <- list.files(paste(tidal.folder, list.dirs.tidal[f], "/", list.dirs.tidal[f], "_result", sep=""), 
					pattern = "*Inserts_Annotated.txt")
	temp <- read.table(file=paste(parent.f, list.dirs.temp[f], "/", input.temp, sep=""), header=T)
	tidal <- read.csv(file=paste(tidal.folder, "/", list.dirs.tidal[f],  "/", list.dirs.tidal[f], 
					"_result", "/", input.tidal, sep=""), 
					sep="\t", header=T)
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
			ind <- which(data1n2[ ,1] == dl[i ,1] & data1n2[ ,2] == dl[i ,2] & data1n2[ ,4] == dl[i ,4] | data1n2[ ,1] == dl[i ,1] & data1n2[ ,3] == dl[i ,3] & data1n2[ ,4] == dl[i ,4])
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

	write.table(mndu2.1, file=paste(parent.r, name, "_new_intervals_1p1.f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
}

### Compare each CONTROL and SELECTION pool to get the frequencies of the intersected intervals and perform Fisher's exact test and Bonferroni correction
for (f in 1:3) {										
	input.c <- read.table(file=paste(parent.r, "CSD",f,"_new_intervals_1p1.f0.10.txt", sep=""), sep="\t", header=T)
	input.r <- read.table(file=paste(parent.r, "RSD",f,"_new_intervals_1p1.f0.10.txt", sep=""), sep="\t", header=T)
	
	input.c.n <- input.c %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "_", remove = TRUE)
	input.r.n <- input.r %>% unite(Chr_TE_ID, c(Chr, TE_ID), sep = "_", remove = TRUE)
	
	data1 <- NULL
	for (i in 1:length(input.c.n[,1])) {
		for (j in 1:length(input.r.n[,1])) {
			if (input.c.n[i,1] == input.r.n[j,1]) {
				if (c(input.c.n[i,2],input.c.n[i,3]) %overlaps% c(input.r.n[j,2], input.r.n[j,3])) {
					pval <- fisher.test(matrix(c(as.numeric(input.c.n[i,5]), as.numeric(input.c.n[i,6]), 
										as.numeric(input.r.n[j,5]), as.numeric(input.r.n[j,6])), nrow=2, 
										ncol=2, byrow=T))$p.val
					data1 <- rbind(data1, c(strsplit(input.c.n[i,1], split="_")[[1]][1], 
									strsplit(input.c.n[i,1], split="_")[[1]][2], input.c.n[i,2], 
									input.c.n[i,3], input.r.n[j,2], input.r.n[j,3], input.c.n[i,4], 
									input.r.n[j,4], scientific(as.numeric(pval), digits=2)))
				}
			}		
		}
	}
	colnames(data1) <- c("Chr", "TE_ID", "Start_C", "Stop_C", "Start_R", "Stop_R", "Freq_C", "Freq_R", "Pval")
	write.table(data1, file=paste(parent.r, "Freqs_CDS-RDS", f,"_f0.10.Pval.txt", sep=""), 
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
		colnames(data2) <- c("Chr", "TE_ID", "Start_C", "Stop_C", "Start_R", "Stop_R", "Freq_C", "Freq_R", "Pval", "Pval_BF")
		write.table(data2, file=paste(parent.r, "Freqs_CDS-RDS", f,"_f0.10.Pval_BF.txt", sep=""), 
					quote=F, sep="\t", row.names=F, col.names=T)
	}
}

### Comparing selection pools
input.1 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS1_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)
input.2 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS2_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)
input.3 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS3_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)

bf.input.1 <- input.1[input.1[, 10] == "**" | input.1[, 10] == "*", ]
bf.input.2 <- input.2[input.2[, 10] == "**" | input.2[, 10] == "*", ]
bf.input.3 <- input.3[input.3[, 10] == "**" | input.3[, 10] == "*", ]

# For each control-selection experiment, get the intervals with TEs corrected by BF with an increase or decrease in frequency from control to selection
bf.input.1.i <- bf.input.1[bf.input.1[,8] > bf.input.1[,7], ]
bf.input.1.d <- bf.input.1[bf.input.1[,8] < bf.input.1[,7], ]
bf.input.2.i <- bf.input.2[bf.input.2[,8] > bf.input.2[,7], ]
bf.input.2.d <- bf.input.2[bf.input.2[,8] < bf.input.2[,7], ]
bf.input.3.i <- bf.input.3[bf.input.3[,8] > bf.input.3[,7], ]
bf.input.3.d <- bf.input.3[bf.input.3[,8] < bf.input.3[,7], ]

write.table(bf.input.1.i, file=paste(parent.r, "Selected_TEs_CDS-RDS1_increase_f_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
write.table(bf.input.1.d, file=paste(parent.r, "Selected_TEs_CDS-RDS1_decrease_f_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
write.table(bf.input.2.i, file=paste(parent.r, "Selected_TEs_CDS-RDS2_increase_f_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
write.table(bf.input.2.d, file=paste(parent.r, "Selected_TEs_CDS-RDS2_decrease_f_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
write.table(bf.input.3.i, file=paste(parent.r, "Selected_TEs_CDS-RDS3_increase_f_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
write.table(bf.input.3.d, file=paste(parent.r, "Selected_TEs_CDS-RDS3_decrease_f_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)

# Get the intervals (Start/Stop) that overlap between the three selected pools for which the frequency changes (between the control and selected) in the same direction (increase or decrease)
data1 <- data2 <- data3 <- NULL
for (i in 1:length(bf.input.1[,1])) {
	for (j in 1:length(bf.input.2[,1])) {
			if (bf.input.1[i,1] == bf.input.2[j,1] && bf.input.1[i,2] == bf.input.2[j,2]) {
				if (c(bf.input.1[i,5],bf.input.1[i,6]) %overlaps% c(bf.input.2[j,5], bf.input.2[j,6])) {
					if (bf.input.1[i,8] > bf.input.1[i,7] && bf.input.2[j,8] > bf.input.2[j,7]) {
						data1 <- bind_cols(bf.input.1[i,1:10], bf.input.2[j,1:10], .id = NULL)
						}
					if (bf.input.1[i,8] < bf.input.1[i,7] && bf.input.2[j,8] < bf.input.2[j,7]) {
						data1 <-  bind_cols(bf.input.1[i,1:10], bf.input.2[j,1:10], .id = NULL)
						}
				}			
			}
			data2 <- rbind(data2, data1)
	}
}
		
data3 <- data2[!duplicated(data2[1:20]), ]

if (!is.null(data3)) {
	data1 <- data2 <- data4 <- NULL
	for (i in 1:length(data3[,1])) {
		for (j in 1:length(bf.input.3[,1])) {
				if (data3[i,11] == bf.input.3[j,1] && data3[i,12] == bf.input.3[j,2]) {
					if (c(data3[i,15],data3[i,16]) %overlaps% c(bf.input.3[j,5], bf.input.3[j,6])) {
						if (data3[i,18] > data3[i,17] && bf.input.3[j,8] > bf.input.3[j,7]) {
							data1 <- bind_cols(data3[i,1:20], bf.input.3[j,1:10], .id = NULL)
							}
						if (data3[i,18] < data3[i,17] && bf.input.3[j,8] < bf.input.3[j,7]) {
							data1 <- bind_cols(data3[i,1:20], bf.input.3[j,1:10], .id = NULL)
							}
					}			
				}
				data2 <- rbind(data2, data1)
		}
	}		
	data4 <- data2[!duplicated(data2[1:30]), ]
	
	if (!is.null(data4)) {
		colnames(data4) <- rep(c("Chr","TE_ID","Start_C","Stop_C","Start_R","Stop_R","Freq_C","Freq_R","Pval","Pval_BF"),3)
		write.table(data4, file=paste(parent.r, "Selected_TEs_all_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
	} else {
		write.table("There are no shared significant TEs between the three pools in which their frequency change in the same direction", 
			file=paste(parent.r, "Selected_TEs_all_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)	
	}
} else {
	write.table("There are no shared significant TEs between the three pools in which their frequency change in the same direction", 
			file=paste(parent.r, "Selected_TEs_all_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}

### Get the significant TEs that increae in frequency only in two selection pools
## Significant TEs after BF correction in RDS1 and RDS2 but not significant in RDS3 (the change in the frequency is in the same direction in the three selective pools)
parent.r <- "/Users/gonzalezlab/Desktop/Analysis/ER_remolina/Results/"

input.1 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS1_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)
input.2 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS2_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)
input.3 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS3_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)

bf.input.1 <- input.1[input.1[, 10] == "**" | input.1[, 10] == "*", ]
bf.input.2 <- input.2[input.2[, 10] == "**" | input.2[, 10] == "*", ]
bf.input.3 <- input.3[input.3[, 10] != "**" & input.3[, 10] != "*", ]

# Get the intervals (Start/Stop) that overlap between the three selected pools for which the frequency changes (between the control and selected) in the same direction (increase or decrease)
data1 <- data2 <- data3 <- NULL
for (i in 1:length(bf.input.1[,1])) {
	for (j in 1:length(bf.input.2[,1])) {
			if (bf.input.1[i,1] == bf.input.2[j,1] && bf.input.1[i,2] == bf.input.2[j,2]) {
				if (c(bf.input.1[i,5],bf.input.1[i,6]) %overlaps% c(bf.input.2[j,5], bf.input.2[j,6])) {
					if (bf.input.1[i,8] > bf.input.1[i,7] && bf.input.2[j,8] > bf.input.2[j,7]) {
						data1 <- bind_cols(bf.input.1[i,1:10], bf.input.2[j,1:10], .id = NULL)
						}
					if (bf.input.1[i,8] < bf.input.1[i,7] && bf.input.2[j,8] < bf.input.2[j,7]) {
						data1 <-  bind_cols(bf.input.1[i,1:10], bf.input.2[j,1:10], .id = NULL)
						}
				}			
			}
			data2 <- rbind(data2, data1)
	}
}		
data3 <- data2[!duplicated(data2[1:20]), ]

if (!is.null(data3)) {
	data1 <- data2 <- data4 <- NULL
	for (i in 1:length(data3[,1])) {
		for (j in 1:length(bf.input.3[,1])) {
				if (data3[i,11] == bf.input.3[j,1] && data3[i,12] == bf.input.3[j,2]) {
					if (c(data3[i,15],data3[i,16]) %overlaps% c(bf.input.3[j,5], bf.input.3[j,6]) && c(data3[i,5],data3[i,6]) %overlaps% c(bf.input.3[j,5], bf.input.3[j,6])) {
						if (data3[i,18] > data3[i,17] && bf.input.3[j,8] > bf.input.3[j,7]) {
							data1 <- bind_cols(data3[i,1:20], bf.input.3[j,1:10], .id = NULL)
							}
						if (data3[i,18] < data3[i,17] && bf.input.3[j,8] < bf.input.3[j,7]) {
							data1 <- bind_cols(data3[i,1:20], bf.input.3[j,1:10], .id = NULL)
							}
					}			
				}
				data2 <- rbind(data2, data1)
		}
	}		
	data4 <- data2[!duplicated(data2[1:30]), ]
	if (!is.null(data4)) {
		colnames(data4) <- rep(c("Chr","TE_ID","Start_C","Stop_C","Start_R","Stop_R","Freq_C","Freq_R","Pval","Pval_BF"),3)
		write.table(data4, file=paste(parent.r, "Selected_TEs_RDS1-RDS2_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
	} else {
		write.table("There are no shared significant TEs between the three pools in which their frequency change in the same direction", 
					file=paste(parent.r, "Selected_TEs_RDS1-RDS2_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
	}	
} else {
	write.table("There are no shared significant TEs between the three pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_RDS1-RDS2_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}

## Significant TEs after BF correction in RDS1 and RDS3 but not significant in RDS2 (the change in the frequency is in the same direction in the three selective pools)
input.1 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS1_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)
input.2 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS3_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)
input.3 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS2_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)

bf.input.1 <- input.1[input.1[, 10] == "**" | input.1[, 10] == "*", ]
bf.input.2 <- input.2[input.2[, 10] == "**" | input.2[, 10] == "*", ]
bf.input.3 <- input.3[input.3[, 10] != "**" & input.3[, 10] != "*", ]

# Get the intervals (Start/Stop) that overlap between the three selected pools for which the frequency changes (between the control and selected) in the same direction (increase or decrease)
data1 <- data2 <- data3 <- NULL
for (i in 1:length(bf.input.1[,1])) {
	for (j in 1:length(bf.input.2[,1])) {
			if (bf.input.1[i,1] == bf.input.2[j,1] && bf.input.1[i,2] == bf.input.2[j,2]) {
				if (c(bf.input.1[i,5],bf.input.1[i,6]) %overlaps% c(bf.input.2[j,5], bf.input.2[j,6])) {
					if (bf.input.1[i,8] > bf.input.1[i,7] && bf.input.2[j,8] > bf.input.2[j,7]) {
						data1 <- bind_cols(bf.input.1[i,1:10], bf.input.2[j,1:10], .id = NULL)
						}
					if (bf.input.1[i,8] < bf.input.1[i,7] && bf.input.2[j,8] < bf.input.2[j,7]) {
						data1 <-  bind_cols(bf.input.1[i,1:10], bf.input.2[j,1:10], .id = NULL)
						}
				}			
			}
			data2 <- rbind(data2, data1)
	}
}		
data3 <- data2[!duplicated(data2[1:20]), ]

if (!is.null(data3)) {
	data1 <- data2 <- data4 <- NULL
	for (i in 1:length(data3[,1])) {
		for (j in 1:length(bf.input.3[,1])) {
				if (data3[i,11] == bf.input.3[j,1] && data3[i,12] == bf.input.3[j,2]) {
					if (c(data3[i,15],data3[i,16]) %overlaps% c(bf.input.3[j,5], bf.input.3[j,6]) && c(data3[i,5],data3[i,6]) %overlaps% c(bf.input.3[j,5], bf.input.3[j,6])) {
						if (data3[i,18] > data3[i,17] && bf.input.3[j,8] > bf.input.3[j,7]) {
							data1 <- bind_cols(data3[i,1:20], bf.input.3[j,1:10], .id = NULL)
							}
						if (data3[i,18] < data3[i,17] && bf.input.3[j,8] < bf.input.3[j,7]) {
							data1 <- bind_cols(data3[i,1:20], bf.input.3[j,1:10], .id = NULL)
							}
					}			
				}
				data2 <- rbind(data2, data1)
		}
	}		
	data4 <- data2[!duplicated(data2[1:30]), ]
	if (!is.null(data4)) {
		colnames(data4) <- rep(c("Chr","TE_ID","Start_C","Stop_C","Start_R","Stop_R","Freq_C","Freq_R","Pval","Pval_BF"),3)
		write.table(data4, file=paste(parent.r, "Selected_TEs_RDS1-RDS3_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
	} else {
	write.table("There are no shared significant TEs between the three pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_RDS1-RDS3_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
	}
} else {
	write.table("There are no shared significant TEs between the three pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_RDS1-RDS3_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}

## Significant TES after BF correction in RDS2 and RDS3 but not significant in RDS1 (the change in the frequency is in the same direction in the three selective pools)
input.1 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS2_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)
input.2 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS3_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)
input.3 <- read.table(file=paste(parent.r, "Freqs_CDS-RDS1_f0.10.Pval_BF.txt", sep=""), sep="\t", header=T)

bf.input.1 <- input.1[input.1[, 10] == "**" | input.1[, 10] == "*", ]
bf.input.2 <- input.2[input.2[, 10] == "**" | input.2[, 10] == "*", ]
bf.input.3 <- input.3[input.3[, 10] != "**" & input.3[, 10] != "*", ]

# Get the intervals (Start/Stop) that overlap between the three selected pools for which the frequency changes (between the control and selected) in the same direction (increase or decrease)
data1 <- data2 <- data3 <- NULL
for (i in 1:length(bf.input.1[,1])) {
	for (j in 1:length(bf.input.2[,1])) {
			if (bf.input.1[i,1] == bf.input.2[j,1] && bf.input.1[i,2] == bf.input.2[j,2]) {
				if (c(bf.input.1[i,5],bf.input.1[i,6]) %overlaps% c(bf.input.2[j,5], bf.input.2[j,6])) {
					if (bf.input.1[i,8] > bf.input.1[i,7] && bf.input.2[j,8] > bf.input.2[j,7]) {
						data1 <- bind_cols(bf.input.1[i,1:10], bf.input.2[j,1:10], .id = NULL)
						}
					if (bf.input.1[i,8] < bf.input.1[i,7] && bf.input.2[j,8] < bf.input.2[j,7]) {
						data1 <-  bind_cols(bf.input.1[i,1:10], bf.input.2[j,1:10], .id = NULL)
						}
				}			
			}
			data2 <- rbind(data2, data1)
	}
}		
data3 <- data2[!duplicated(data2[1:20]), ]

if (!is.null(data3)) {
	data1 <- data2 <- data4 <- NULL
	for (i in 1:length(data3[,1])) {
		for (j in 1:length(bf.input.3[,1])) {
				if (data3[i,11] == bf.input.3[j,1] && data3[i,12] == bf.input.3[j,2]) {
					if (c(data3[i,15],data3[i,16]) %overlaps% c(bf.input.3[j,5], bf.input.3[j,6]) && c(data3[i,5],data3[i,6]) %overlaps% c(bf.input.3[j,5], bf.input.3[j,6])) {
						if (data3[i,18] > data3[i,17] && bf.input.3[j,8] > bf.input.3[j,7]) {
							data1 <- bind_cols(data3[i,1:20], bf.input.3[j,1:10], .id = NULL)
							}
						if (data3[i,18] < data3[i,17] && bf.input.3[j,8] < bf.input.3[j,7]) {
							data1 <- bind_cols(data3[i,1:20], bf.input.3[j,1:10], .id = NULL)
							}
					}			
				}
				data2 <- rbind(data2, data1)
		}
	} 		
	data4 <- data2[!duplicated(data2[1:30]), ]
	if (!is.null(data4)) {
		colnames(data4) <- rep(c("Chr","TE_ID","Start_C","Stop_C","Start_R","Stop_R","Freq_C","Freq_R","Pval","Pval_BF"),3)
		write.table(data4, file=paste(parent.r, "Selected_TEs_RDS2-RDS3_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
	} else {
		write.table("There are no shared significant TEs between the three pools in which their frequency change in the same direction",
					file=paste(parent.r, "Selected_TEs_RDS2-RDS3_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
	}
} else {
	write.table("There are no shared significant TEs between the three pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_RDS2-RDS3_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}

### Private significant TEs
input.1 <- read.table(file=paste(parent.r, "Selected_TEs_CDS-RDS1_increase_f_BF_f0.10.txt", sep=""),
						sep="\t", header=T)
input.2 <- read.table(file=paste(parent.r, "Selected_TEs_CDS-RDS2_increase_f_BF_f0.10.txt", sep=""),
						sep="\t", header=T)
input.3 <- read.table(file=paste(parent.r, "Selected_TEs_CDS-RDS3_increase_f_BF_f0.10.txt", sep=""),
						sep="\t", header=T)

control.1 <- read.table(file=paste(parent.r, "CSD1_new_intervals_1p1.f0.10.txt", sep=""), sep="\t", header=T)
control.2 <- read.table(file=paste(parent.r, "CSD2_new_intervals_1p1.f0.10.txt", sep=""), sep="\t", header=T) 
control.3 <- read.table(file=paste(parent.r, "CSD3_new_intervals_1p1.f0.10.txt", sep=""), sep="\t", header=T) 

control.12 <- rbind(control.1, control.2)
control.13 <- rbind(control.1, control.3)
control.23 <- rbind(control.2, control.3)

control.tm1 <- read.table(file=paste(parent.f, "CSD1-deSim-mel/CSD1.insertion.1p1.final.txt", sep=""), sep="\t", header=T)
control.tm2 <- read.table(file=paste(parent.f, "CSD2-deSim-mel/CSD2.insertion.1p1.final.txt", sep=""), sep="\t", header=T) 
control.tm3 <- read.table(file=paste(parent.f, "CSD3-deSim-mel/CSD3.insertion.1p1.final.txt", sep=""), sep="\t", header=T)

control.td1 <- read.csv(file=paste(tidal.folder, "CSD1/insertion/CSD1_Inserts_Annotated.txt", sep=""), sep="\t", header=T)
control.td2 <- read.csv(file=paste(tidal.folder, "CSD2/insertion/CSD2_Inserts_Annotated.txt", sep=""), sep="\t", header=T) 
control.td3 <- read.csv(file=paste(tidal.folder, "CSD3/insertion/CSD3_Inserts_Annotated.txt", sep=""), sep="\t", header=T)

control.td1sub <- control.td1[,2:5]
control.td2sub <- control.td2[,2:5]
control.td3sub <- control.td3[,2:5]

colnames(control.td1sub) <- c("Chr","Start", "End",	"TransposonName")
colnames(control.td2sub) <- c("Chr","Start", "End",	"TransposonName")
colnames(control.td3sub) <- c("Chr","Start", "End",	"TransposonName")

control.tt.12 <- rbind(control.tm1[,1:4], control.td1sub, control.tm2[,1:4], control.td2sub)
control.tt.13 <- rbind(control.tm1[,1:4], control.td1sub, control.tm3[,1:4], control.td3sub)
control.tt.23 <- rbind(control.tm2[,1:4], control.td2sub, control.tm3[,1:4], control.td3sub)

# Significant TEs in CDS-RDS1 but not present in CSD2 and CSD3 using the intevals from the TEMP and TIDAL output files
data1 <- data2 <- NULL
for (i in 1:length(input.1[,1])) {
	data1 <- NULL
	for (j in 1:length(control.tt.23[,1])) {
			if (input.1[i,1] == control.tt.23[j,1] && input.1[i,2] == control.tt.23[j,4]) {
				if (c(input.1[i,3],input.1[i,4]) %overlaps% c(control.tt.23[j,2], control.tt.23[j,3])) {
					data1 <- rbind(data1, input.1[i,])
				}	
			}		
	}	
	if (is.null(data1)) {
		data2 <- rbind(data2, input.1[i,])
	}
}

if (!is.null(data2)) {
	write.table(data2, file=paste(parent.r, "Selected_TEs_CDS-RDS1_unique_TT_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
} else {	
	write.table("There are no shared TEs between pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_CDS-RDS1_unique_TT_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
}

# Significant TEs in CDS-RDS1 but not present in CSD2 and CSD3 using the new intevals
data1 <- data2 <- NULL
for (i in 1:length(input.1[,1])) {
	data1 <- NULL
	for (j in 1:length(control.23[,1])) {
			if (input.1[i,1] == control.23[j,1] && input.1[i,2] == control.23[j,4]) {
				if (c(input.1[i,3],input.1[i,4]) %overlaps% c(control.23[j,2], control.23[j,3])) {
					data1 <- rbind(data1, input.1[i,])
				}	
			}		
	}	
	if (is.null(data1)) {
		data2 <- rbind(data2, input.1[i,])
	}
}	

if (!is.null(data2)) {
	write.table(data2, file=paste(parent.r, "Selected_TEs_CDS-RDS1_unique_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
} else {	
	write.table("There are no shared TEs between pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_CDS-RDS1_unique_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
}

# Significant TEs in CDS-RDS2 against CSD1 and CSD3 using the intevals from the TEMP and TIDAL output files
data1 <- data2 <- NULL
for (i in 1:length(input.2[,1])) {
	data1 <- NULL
	for (j in 1:length(control.tt.13[,1])) {
			if (input.2[i,1] == control.tt.13[j,1] && input.2[i,2] == control.tt.13[j,4]) {
				if (c(input.2[i,3],input.2[i,4]) %overlaps% c(control.tt.13[j,2], control.tt.13[j,3])) {
					data1 <- rbind(data1, input.2[i,])
				}	
			}		
	}	
	if (is.null(data1)) {
		data2 <- rbind(data2, input.2[i,])
	}
}	

if (!is.null(data2)) {
	write.table(data2, file=paste(parent.r, "Selected_TEs_CDS-RDS2_unique_TT_BF_f0.10.txt", sep=""),
				quote=F, sep="\t", row.names=F, col.names=T)
} else {	
	write.table("There are no significant TEs between pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_CDS-RDS2_unique_TT_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
}

# Significant TEs in CDS-RDS2 against CSD1 and CSD3 using the new intevals
data1 <- data2 <- NULL
for (i in 1:length(input.2[,1])) {
	data1 <- NULL
	for (j in 1:length(control.13[,1])) {
			if (input.2[i,1] == control.13[j,1] && input.2[i,2] == control.13[j,4]) {
				if (c(input.2[i,3],input.2[i,4]) %overlaps% c(control.13[j,2], control.13[j,3])) {
					data1 <- rbind(data1, input.2[i,])
				}	
			}		
	}	
	if (is.null(data1)) {
		data2 <- rbind(data2, input.2[i,])
	}
}	

if (!is.null(data2)) {
	write.table(data2, file=paste(parent.r, "Selected_TEs_CDS-RDS2_unique_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
} else {	
	write.table("There are no shared TEs between pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_CDS-RDS2_unique_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
}

# Significant TEs in CDS-RDS3 against CSD1 and CSD2 using the intevals from the TEMP and TIDAL output files
data1 <- data2 <- NULL
for (i in 1:length(input.3[,1])) {
	data1 <- NULL
	for (j in 1:length(control.tt.12[,1])) {
			if (input.3[i,1] == control.tt.12[j,1] && input.3[i,2] == control.tt.12[j,4]) {
				if (c(input.3[i,3],input.3[i,4]) %overlaps% c(control.tt.12[j,2], control.tt.12[j,3])) {
					data1 <- rbind(data1, input.3[i,])
				}	
			}		
	}	
	if (is.null(data1)) {
		data2 <- rbind(data2, input.3[i,])
	}
}	

if (!is.null(data2)) {
	write.table(data2, file=paste(parent.r, "Selected_TEs_CDS-RDS3_unique_TT_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
} else {	
	write.table("There are no sharedTEs between pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_CDS-RDS3_unique_TT_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
}

# Significant TEs in CDS-RDS3 against CSD1 and CSD2 using new intervals
data1 <- data2 <- NULL
for (i in 1:length(input.3[,1])) {
	data1 <- NULL
	for (j in 1:length(control.12[,1])) {
			if (input.3[i,1] == control.12[j,1] && input.3[i,2] == control.12[j,4]) {
				if (c(input.3[i,3],input.3[i,4]) %overlaps% c(control.12[j,2], control.12[j,3])) {
					data1 <- rbind(data1, input.3[i,])
				}	
			}		
	}	
	if (is.null(data1)) {
		data2 <- rbind(data2, input.3[i,])
	}
}	

if (!is.null(data2)) {
	write.table(data2, file=paste(parent.r, "Selected_TEs_CDS-RDS3_unique_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
} else {	
	write.table("There are no shared TEs between pools in which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_CDS-RDS3_unique_BF_f0.10.txt", sep=""), 
				quote=F, sep="\t", row.names=F, col.names=T)
}
