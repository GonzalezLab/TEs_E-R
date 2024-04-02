#!/usr/bin/env Rscript

library(GSEABase)
library(GOstats)
library("org.Dm.eg.db")
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
#args[1]: path to the folder with the diferent subfolders (pools) containing the TEMP output files (i.e., "/home/")
#args[2]: path to the results folder
#args[3]: name of the file containing the candidate genes (columns are: Gene_name; FlyBase_gene_ID)
#args[4]: name of the output file without extension (e.g, "results_remolina_go_analysis")  							

parent.f <- args[1]
parent.o <- args[2]
name <- args[3]
out <- args[4]

# Create a goAllFrame object
frame = toTable(org.Dm.egGO)
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
goFrame = GOFrame(goframeData, organism="Drosophila melanogaster")
goAllFrame=GOAllFrame(goFrame)

gsc <- GeneSetCollection(goAllFrame, setType=GOCollection())

# Get the genes in the universe
universe = Lkeys(org.Dm.egGO)
print("Universer load")

# Read the file with the candidate genes
candidates <- read.csv(file=paste(parent.f, name, sep=""), header=T, sep="\t")

# Get the entrex ID of the candidates genes from their Flybase IDs
a <- org.Dm.egFLYBASE
mapped_genes <- mappedkeys(a)
b <- as.data.frame(a[mapped_genes])

genes <- b[b[,2] %in% candidates[, 2], 1]

# Set a specific P-vlaue cutoff
pval <- 0.01

# Set the parameters
params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                            geneSetCollection = gsc,
                            geneIds = genes,
                            universeGeneIds = universe,
                            ontology = "BP",
                            pvalueCutoff = pval,
                            conditional = TRUE,
                            testDirection = "over")

# Perform GO enrichment (Biological Process)                           
Over_BP <- hyperGTest(params)
write.table(summary(Over_BP), file = paste(parent.o, out, "_BP.txt", sep=""), 
			quote=FALSE, col.names=T, row.names=F, sep="\t")

# Perform GO enrichment (Molecular Function)
ontology(params) <- "MF"
Over_MF <- hyperGTest(params)	
write.table(summary(Over_MF), file = paste(parent.o, out, "_MF.txt", sep=""), 
			quote=FALSE, col.names=T, row.names=F, sep="\t")

# # Perform GO enrichment (Cellular Component)
ontology(params) <- "CC"
Over_CC <- hyperGTest(params)	
write.table(summary(Over_CC), file = paste(parent.o, out, "_CC.txt", sep=""), 
			quote=FALSE, col.names=T, row.names=F, sep="\t")		

## Correct for multiple testing
# BP
params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                            geneSetCollection = gsc,
                            geneIds = genes,
                            universeGeneIds = universe,
                            ontology = "BP",
                            pvalueCutoff = 1,
                            conditional = TRUE,
                            testDirection = "over")
                            
Over_BP <- hyperGTest(params)
over_BPqval <- qvalue(summary(Over_BP)$Pvalue, fdr.level=0.1, lambda=0, pi0.method="smoother")

# Write the q-values results
write.qvalue(over_BPqval, file = paste(parent.o, out, "_BP_qval.txt", sep=""), 
			sep = " ", eol = "\n", na = "NA",
            row.names = FALSE, col.names = TRUE)
# MF			
ontology(params) <- "MF"
Over_MF <- hyperGTest(params)
over_MFqval <- qvalue(summary(Over_MF)$Pvalue, fdr.level=0.1, lambda=0, pi0.method="smoother")

# Write the q-values results
write.qvalue(over_MFqval, file = paste(parent.o, out, "_MF_qval.txt", sep=""), 
			sep = " ", eol = "\n", na = "NA",
            row.names = FALSE, col.names = TRUE)
			
# CC			
ontology(params) <- "CC"
Over_CC <- hyperGTest(params)
over_CCqval <- qvalue(summary(Over_CC)$Pvalue, fdr.level=0.1, lambda=0, pi0.method="smoother")

# Write the q-values results
write.qvalue(over_CCqval, file = paste(parent.o, out, "_CC_qval.txt", sep=""), 
			sep = " ", eol = "\n", na = "NA",
            row.names = FALSE, col.names = TRUE)		
            