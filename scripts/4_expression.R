#
# Call expression from aligned SRA reads
# copyright (c) 2022 - Danny Arends
#

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")
library("vioplot")
library("RColorBrewer")
library("biomaRt")

# Go into the output folder
setwd("/deac/inf/adminGrp/anderss/scratch/zhangGrp")

# Create DB and exons per Gene
db <- makeTxDbFromGFF("genome/Saccharomyces_cerevisiae.R64-1-1.108.gtf", 
                      format = "gtf", organism = "Saccharomyces", 
                      dataSource = "https://ftp.ensembl.org/pub/release-108/")

# Get the exons per gene, and compute bp lengths of all genes
exons <- exonsBy(db, by = "gene")
gene.lengths <- lapply(exons, function(x){ sum(width(reduce(x))) })

setwd("/deac/inf/adminGrp/anderss/scratch/zhangGrp/output")

# Samples and simple names
samples <- c("SRR13978640", "SRR13978641", "SRR13978642", "SRR13978643", "SRR13978644", "SRR13978645")
names(samples) <- c("SPRC_1", "SPRC_2", "SPRC_3", "CTRL_1", "CTRL_2", "CTRL_3")

# Check if BAM files for all samples exist
files <- c()
for(s in samples){
  fp <- paste0(s, ".aln/", s, "Aligned.sortedByCoord.RD.RG.RC.out.bam")
  if(file.exists(fp)) files <- c(files, fp)
}

# Load in the BAM files
bams <- BamFileList(files, yieldSize = 100000, asMates=TRUE)

# Overlap BAM reads and genes
overlap <- summarizeOverlaps(exons, bams, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

# Extract the raw-reads per gene
readcount <- assay(overlap)
colnames(readcount) <- gsub("Aligned.sortedByCoord.RD.RG.RC.out.bam", "", colnames(readcount))
write.table(readcount, "readcount.raw.txt", sep = "\t", quote = FALSE)

# Calculate the RPKM values per gene
# RPKM = (10^9 * C)/(N * L)
# C = Number of reads mapped to a gene
# N = Total mapped reads in the sample
# L = gene length in base-pairs for a gene

# Get the total number of reads per samples
N <- apply(readcount, 2, sum)

# Loop through all genes, compute RPKM
n <- 1
RPKM <- t(apply(readcount, 1, function(C){
  L     <- as.numeric(gene.lengths[n])
  RPKM  <- (10^9 * C) / (N * L)
  n    <<- n + 1
  return(round(RPKM, d = 1))
}))

op <- par(mar = c(8,4,2,2))
#Violin distribution plot
vioplot(RPKM, col = c(rep("orange", 3), rep("lightblue", 3)), ylab = "reads", las = 2)
legend("topleft", c("SPRC", "CTRL"), fill = c("orange", "lightblue"))

write.table(RPKM, "RPKM.txt", sep = "\t", quote = FALSE)

q("no")

# Quantile normalization of RPKM values
RPKM.norm <- round(normalize.quantiles(as.matrix(RPKM)), d = 1)
colnames(RPKM.norm) <- colnames(RPKM)
rownames(RPKM.norm) <- rownames(RPKM)

#Violin distribution plot
vioplot(RPKM.norm, col = c(rep("orange", 3), rep("lightblue",3)), ylab = "normalized", las = 2)
legend("topleft", c("SPRC", "CTRL"), fill = c("orange", "lightblue"))

write.table(RPKM.norm, "RPKM.norm.txt", sep = "\t", quote = FALSE)

# LOG2 Qnorm RPKM (so we can treat it as microarray data)
RPKM.l2 <- round(log2(RPKM.norm), d = 1)
RPKM.l2[RPKM.l2 < 0] <- 0

#Violin distribution plot
vioplot(RPKM.l2, col = c(rep("orange", 3), rep("lightblue",3)), ylab = "log2(normalize)", las = 2)
legend("topleft", c("SPRC", "CTRL"), fill = c("orange", "lightblue"))

write.table(RPKM.l2, "RPKM.norm.log2.txt", sep = "\t", quote = FALSE)

# P-values and Log2 fold change
pvals <- apply(RPKM.l2, 1, function(x){
  tryCatch(t.test(x[1:3], x[4:6])$p.value, error = function(x){return(NA);})
})

fc <- apply(RPKM.l2, 1, function(x){
  tryCatch(log2(mean(x[1:3]) / mean(x[4:6])), error = function(x){return(NA);})
})

# Assign colors based on P-values
colz <- rep("black", length(pvals))
colz[which(pvals < 5e-2)] <- "red"
colz[which(pvals < 1e-2)] <- "gold"
colz[which(pvals < 1e-3)] <- "blue"

# Volcano plot (x = fc, y = -log10(P-values))
plot(fc, log10(pvals), col=colz, pch=18, main ="Vulcano plot", xlab="Fold Change")
legend("topleft", pch=18, c("<0.05", "<0.01", "<0.001"), col = c("red", "gold", "blue"))

# Down & Up regulated genes
down <- RPKM.l2[which(pvals < 5e-2 & fc < -0.3),]
dclust <- down[hclust(dist(down))$order,]
up <- RPKM.l2[which(pvals < 5e-2 & fc > 0.3),]
uclust <- up[hclust(dist(up))$order,]

# Gene IDs of up/Down regulated genes
geneIDs <- c(rownames(dclust), rownames(uclust))

# Custom heatmap using the spectral colors
op <- par(mar = c(8, 6, 2,1))
image(x = 1:ncol(RPKM.l2), 
      y = 1:(nrow(down)+nrow(up)), 
      z = t(rbind(dclust, uclust)), 
      xaxt='n',yaxt='n',xlab="", ylab="", 
      col = brewer.pal(11, "Spectral"))

axis(2, at = 1:length(geneIDs), labels = geneIDs, las = 2, cex.axis=0.7)
axis(1, at = 1:3, labels = colnames(RPKM.l2)[1:3], las = 2, col.axis = "orange")
axis(1, at = 4:6, labels = colnames(RPKM.l2)[4:6], las = 2, col.axis = "blue")

# Compute mean expression and standard deviations
means <- t(apply(RPKM.l2, 1, function(x){
  tryCatch(round(c(mean(x[1:3]),mean(x[4:6])),1), error = function(x){return(NA);})
}))

sds <- t(apply(RPKM.l2, 1, function(x){
  tryCatch(round(c(sd(x[1:3]),sd(x[4:6])),1), error = function(x){return(NA);})
}))

colnames(means) <- c("SPRC", "CTRL")
colnames(sds) <- c("SPRC", "CTRL")

# Create an overview table
overview <- cbind("CTRL" = means[, "CTRL"], 
                  "CTRL(SD)" = sds[, "CTRL"], 
                  "SPRC" = means[, "SPRC"], 
                  "SPRC(SD)" = sds[, "SPRC"], 
                  FC = round(fc,1),
                  P = round(pvals,6))

overview[1:10,]

# Use biomaRt to retrieve gene names, location, and description
library(biomaRt)
bio.mart <- useMart("ensembl", "scerevisiae_gene_ensembl")



mattr <- c("ensembl_gene_id", "external_gene_name", 
           "chromosome_name", "start_position", "end_position",
           "description")

res.bm <- getBM(attributes = mattr, 
                filters = c("ensembl_gene_id"), 
                values = geneIDs, mart = bio.mart)
rownames(res.bm) <- res.bm[, "ensembl_gene_id"]

# Merge biomaRt results with the overview
p1 <- res.bm[geneIDs, c("external_gene_name", "chromosome_name", "start_position", "end_position")]
overview <- cbind(p1, overview[geneIDs,], res.bm[geneIDs, "description"])
colnames(overview)[1:4] <- c("GeneName", "Chr", "Start", "End")
colnames(overview)[11] <- c("Description")

overview[1:10,1:10]

# Write out the table
write.table(overview, "overview.ann.txt", sep = "\t", quote = FALSE)
