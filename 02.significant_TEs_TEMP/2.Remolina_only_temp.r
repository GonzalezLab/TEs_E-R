#!/usr/bin/env Rscript

library(DescTools) 
library(GenomicRanges)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with the diferent subfolders (pools) containing the TEMP output files (i.e., "/home/")
#args[2]: path to the results folder

### Read TEMP output files (avoid this step if the script 1.Remolina.r was already run before)
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

### Compare each CONTROL and SELECTION pool to get the frequencies of the intersected intervals and perform Fisher's exact test and Bonferroni correction
parent.r <- args[2]

list.dirs <- dir(parent.f, pattern ="*-deSim-mel")

for (d in 1:length(list.dirs.temp)) {
	temp.f <- NULL
	input.temp <- list.files(paste(parent.f, list.dirs[d], sep=""), pattern = "*insertion.1p1.f0.10.txt")
	temp <- read.table(file=paste(parent.f, list.dirs[d], "/", input.temp, sep=""), header=T)
	temp.f <- cbind(temp, matrix("NA", ncol=1, nrow=length(temp[,1])))
	temp.f[, 14] <- as.numeric(temp.f[, 14])
	name <- strsplit(input.temp, ".", fixed=T)[[1]][1]
	for (l in 1:length(temp[,1])) {
		ref <- round((as.numeric(temp[l,7]) - (as.numeric(temp[l,8])*as.numeric(temp[l,7])))/as.numeric(temp[l,8]))	
		temp.f[l, 14] <- round((as.numeric(temp[l,7]) - (as.numeric(temp[l,8])*as.numeric(temp[l,7])))/as.numeric(temp[l,8]))	
	}
	colnames(temp.f)[14]<- "Reads_ref"
	write.table(temp.f, file=paste(parent.r, name, ".insertion.1p1.f0.10_ref.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}

for (f in 1:3) {
	input.c <- read.table(file=paste(parent.r, "CSD", f, ".insertion.1p1.f0.10_ref.txt", sep=""), sep="\t", header=T)
	input.r <- read.table(file=paste(parent.r, "RSD", f, ".insertion.1p1.f0.10_ref.txt", sep=""), sep="\t", header=T)
	
	data1 <- NULL
	for (i in 1:length(input.c[,1])) {
		for (j in 1:length(input.r[,1])) {
			if (input.c[i,1] == input.r[j,1] && input.c[i,4] == input.r[j,4]) {
				if (c(input.c[i,2], input.c[i,3]) %overlaps% c(input.r[j,2], input.r[j,3])) {
					pval <- fisher.test(matrix(c(as.numeric(input.c[i,7]), as.numeric(input.c[i,14]), as.numeric(input.r[j,7]), 
										as.numeric(input.r[j,14])), nrow=2, ncol=2, byrow=T))$p.val
					data1 <- rbind(data1, c(as.character(input.c[i,1]), as.character(input.c[i,4]), as.numeric(input.c[i,2]), 
									as.numeric(input.c[i,3]), as.numeric(input.r[j,2]), as.numeric(input.r[j,3]), as.numeric(input.c[i,8]), 
									as.numeric(input.r[j,8]), scientific(as.numeric(pval), digits=2)))
				}
			}		
		}
	}
	colnames(data1) <- c("Chr", "TE_ID", "Start_C", "Stop_C", "Start_R", "Stop_R", "Freq_C", "Freq_R", "Pval")
	write.table(data1, file=paste(parent.r, "Freqs_CDS-RDS", f,"_f0.10.Pval.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)

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
		write.table(data2, file=paste(parent.r, "Freqs_CDS-RDS", f,"_f0.10.Pval_BF.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
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
		write.table(data4, file=paste(parent.r, "Selected_TEs_all_BF_f0.10.txt", sep=""), 
					quote=F, sep="\t", row.names=F, col.names=T)
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
		write.table(data4, file=paste(parent.r, "Selected_TEs_RDS1-RDS2_BF_f0.10txt", sep=""), quote=F, sep="\t", col.names=T, row.names=F)
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
		write.table(data4, file=paste(parent.r, "Selected_TEs_RDS1-RDS3_BF_f0.10.txt", sep=""), 
					quote=F, sep="\t", row.names=F, col.names=T)
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

#This is to extrcat the intervals (start Stop) that overlap between the three selected replicates for which the frequency changes (between the control and selected) in the same direction (inv¡crease or decrease)
#The loops are separated to increase the speed of the process

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
		write.table(data4, file=paste(parent.r, "Selected_TEs_RDS2-RDS3_BF_f0.10txt", sep=""), 
					quote=F, sep="\t", col.names=T, row.names=F)
	} else {
	write.table("There are no shared significant TEs between the three pools which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_RDS2-RDS3_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
	}	
} else {
	write.table("There are no shared significant TEs between the three pools which their frequency change in the same direction", 
				file=paste(parent.r, "Selected_TEs_RDS2-RDS3_BF_f0.10.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}

### Private significant TEs
input.1 <- read.table(file=paste(parent.r, "Selected_TEs_CDS-RDS1_increase_f_BF_f0.10.txt", sep=""), sep="\t", header=T)
input.2 <- read.table(file=paste(parent.r, "Selected_TEs_CDS-RDS2_increase_f_BF_f0.10.txt", sep=""), sep="\t", header=T)
input.3 <- read.table(file=paste(parent.r, "Selected_TEs_CDS-RDS3_increase_f_BF_f0.10.txt", sep=""), sep="\t", header=T)

control.1 <- read.table(file=paste(parent.r, "CSD1.insertion.1p1_ref.txt", sep=""), sep="\t", header=T)
control.2 <- read.table(file=paste(parent.r, "CSD2.insertion.1p1_ref.txt", sep=""), sep="\t", header=T) 
control.3 <- read.table(file=paste(parent.r, "CSD3.insertion.1p1_ref.txt", sep=""), sep="\t", header=T) 

control.12 <- rbind(control.1, control.2)
control.13 <- rbind(control.1, control.3)
control.23 <- rbind(control.2, control.3)

# Significant TEs in CDS-RDS1 but not present in CSD2 and CSD3 using the intevals from the TEMP output files
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

# Significant TEs in CDS-RDS2 but not present in CSD1 and CSD3 using the intevals from the TEMP output files
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

# Significant TEs in CDS-RDS3 but not present in CSD1 and CSD2 using the intevals from the TEMP output files
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