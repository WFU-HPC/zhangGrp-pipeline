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


library(edgeR)
library(xlsx)

###Functions
#read in all genes featureCounts table
readit = function(filename,dfout){
  dfout = read.table(file=filename,header=T,row.names=1,sep="\t")
  names(dfout) = sub("Aligned.sortedByCoord.out.bam","",names(dfout))
  dfout = dfout[,!grepl("^Geneid",names(dfout))]
  dfout = floor(dfout)
  return(dfout)
}

#read in forward featureCounts table
readitfwd = function(filename,dfout){
  dfout = read.table(file=filename,header=T,row.names=1,sep="\t")
  names(dfout) = sub("_forward.bam","",names(dfout))
  dfout = dfout[,!grepl("^Geneid",names(dfout))]
  dfout = floor(dfout)
  return(dfout)
}

#read in reverse featureCounts table
readitrev = function(filename,dfout){
  dfout = read.table(file=filename,header=T,row.names=1,sep="\t")
  names(dfout) = sub("_reverse.bam","",names(dfout))
  dfout = dfout[,!grepl("^Geneid",names(dfout))]
  dfout = floor(dfout)
  return(dfout)
}

#DGE analysis using edgeR
dgplus = function(df,group1,group2,grouping){
  grp1 = group1
  grp2 = group2
  #subset original dataframe by merging string groups
  temp1 = subset(df,select=grp1)
  temp2 = subset(df,select=grp2)
  #new subset dataframe below
  df = transform(merge(temp1,temp2,by="row.names"),row.names=Row.names,Row.names=NULL)
  
  y <- df
  y_full <- DGEList(counts=y,group = grouping)
  keep <- rowSums(cpm(y_full) > 0) >= 2
  lm1.y <- y_full[keep, keep.lib.sizes =FALSE]
  lm1.y <- calcNormFactors(lm1.y, method="TMM")
  lm1.tmm = as.data.frame(lm1.y$counts)
  lm1.design <- model.matrix( ~ as.numeric(group==1), data = lm1.y$samples)
  lm1.v <- voom(lm1.y, lm1.design, plot=FALSE)
  lm1.fit <- lmFit(lm1.v, lm1.design)
  lm1.fit.c <- eBayes(lm1.fit)
  lm1.top <- topTable(lm1.fit.c, sort="P",number='all')
  lm1.top.merge = transform(merge(lm1.top,df,by="row.names"),row.names=Row.names,Row.names=NULL)
  lm1.top.merge = lm1.top.merge[order(lm1.top.merge$adj.P.Val),]
  return(lm1.top.merge)
}

#phased out this function
list2excelout = function(datalist,sheetnamelist,filename){
  wb <- createWorkbook()
  datas <- datalist
  sheetnames = names(sheetnamelist)
  sheets <- lapply(sheetnames, createSheet, wb = wb)
  void <- Map(addDataFrame, datas, sheets)
  saveWorkbook(wb, file = filename)
}

#subsetting significant genes into gene lists for GO
#these files can be input into Ontologizer (Java)
#Right now significance is hard-coded to specific P-Value threshold
list2subsettable = function(datalist,sheetnamelist,fileextension){
  #first element of datalist and sheetname list are empty, that's why -1
  sublist = lapply(datalist[-1], function(x) {subset(x,adj.P.Val < .0001)} )
  for (strain in seq_along(sublist)){
    filename = paste0(names(sheetnamelist[-1])[strain], fileextension)
    write.table( rownames(sublist[[strain]]) ,filename,quote=F,sep="\t",row.names=F,col.names=F)
  }
  
}

#Output the population list of genes for each contrast tested
list2poptable = function(datalist,sheetnamelist,fileextension){
  #first element of datalist and sheetname list are empty, that's why -1
  #sublist = lapply(datalist[-1], function(x) {subset(x,adj.P.Val < .0001)} )
  sublist=datalist[-1]
  for (strain in seq_along(sublist)){
    filename = paste0(names(sheetnamelist[-1])[strain], fileextension)
    write.table( rownames(sublist[[strain]]) ,filename,quote=F,sep="\t",row.names=F,col.names=F)
  }
  
}

#Gene tables
pombe_geneOfrac_all = readit("Spombe_gene_withO_fraction_all.txt")
pombe_geneOfrac_fwd = readitfwd("Spombe_gene_withO_fraction_fwd.txt")
pombe_geneOfrac_rev = readitrev("Spombe_gene_withO_fraction_rev.txt")

nontemp_lst = list(
  "wt" = "A1_WT_1|B1_WT_2",
  "lsd1" = "A3|B3",
  "lsd2" = "A5|B5",
  "clr4" = "A7|B7",
  "clr4_lsd1" = "A8|B8",
  "set1" = "A9|B9",
  "set1_lsd2" = "A11|B11",
  "clr6" = "A12|B12",
  "clr6_lsd1" = "A13|B13",
  "clr6_lsd2" = "A14|B14",
  "sir2" = "A15|B15",
  "sir2_lsd1" = "A16|B16",
  "sir2_lsd2" = "A17|B17"
  
)

temp_lst = list(
  "wt_37" = "A2|B2",
  "lsd1_37" = "A4|B4",
  "lsd2_37" = "A6|B6",
  "set1_37" = "A10|B10"
)

#all genes
nontempdgeout = vector("list",(length(nontemp_lst)-1) )
for (strain in (seq(2:(length(nontemp_lst)) )+1) ){
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_all),value=T)
  mutants = grep(nontemp_lst[[strain]],names(pombe_geneOfrac_all),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  nontempdgeout[[strain]]=dgplus(pombe_geneOfrac_all,controls,mutants,compare_grp)
  
}

#forward strand
nontempdgeout_fwd = vector("list",(length(nontemp_lst)-1) )
for (strain in (seq(2:(length(nontemp_lst)) )+1) ){
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_fwd),value=T)
  mutants = grep(nontemp_lst[[strain]],names(pombe_geneOfrac_fwd),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  nontempdgeout_fwd[[strain]]=dgplus(pombe_geneOfrac_fwd,controls,mutants,compare_grp)
  
}

#reverse strand
nontempdgeout_rev = vector("list",(length(nontemp_lst)-1) )
for (strain in (seq(2:(length(nontemp_lst)) )+1) ){
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_rev),value=T)
  mutants = grep(nontemp_lst[[strain]],names(pombe_geneOfrac_rev),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  nontempdgeout_rev[[strain]]=dgplus(pombe_geneOfrac_rev,controls,mutants,compare_grp)
  
}

#Export into Excel files
wb <- createWorkbook()
datas <- nontempdgeout
sheetnames = names(nontemp_lst)
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, datas, sheets)
saveWorkbook(wb, file = "nontempall_rnaseq.xlsx")

#list2excelout(nontempdgeout_fwd,nontemp_lst,"nontempfwd_rnaseq.xlsx")  
#list2excelout(nontempdgeout_rev,nontemp_lst,"nontemprev_rnaseq.xlsx")  

list2subsettable(nontempdgeout,nontemp_lst,"_nontemp_all_list.txt")
list2subsettable(nontempdgeout_fwd,nontemp_lst,"_nontemp_fwd_list.txt")
list2subsettable(nontempdgeout_rev,nontemp_lst,"_nontemp_rev_list.txt")

#temperature mutants
tempdgeout = vector("list",(length(temp_lst)-1) )
for (strain in (seq(2:(length(temp_lst)) )+1) ){
  controls = grep(temp_lst[[1]],names(pombe_geneOfrac_all),value=T)
  mutants = grep(temp_lst[[strain]],names(pombe_geneOfrac_all),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  tempdgeout[[strain]]=dgplus(pombe_geneOfrac_all,controls,mutants,compare_grp)
  
}

#forward strand
tempdgeout_fwd = vector("list",(length(temp_lst)-1) )
for (strain in (seq(2:(length(temp_lst)) )+1) ){
  controls = grep(temp_lst[[1]],names(pombe_geneOfrac_fwd),value=T)
  mutants = grep(temp_lst[[strain]],names(pombe_geneOfrac_fwd),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  tempdgeout_fwd[[strain]]=dgplus(pombe_geneOfrac_fwd,controls,mutants,compare_grp)
  
}

#reverse strand
tempdgeout_rev = vector("list",(length(temp_lst)-1) )
for (strain in (seq(2:(length(temp_lst)) )+1) ){
  controls = grep(temp_lst[[1]],names(pombe_geneOfrac_rev),value=T)
  mutants = grep(temp_lst[[strain]],names(pombe_geneOfrac_rev),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  tempdgeout_rev[[strain]]=dgplus(pombe_geneOfrac_rev,controls,mutants,compare_grp)
  
}

#output here
wb <- createWorkbook()
datas <- tempdgeout
sheetnames = names(temp_lst)
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, datas, sheets)
saveWorkbook(wb, file = "temp_rnaseq.xlsx")

#list2excelout(tempdgeout_fwd,temp_lst,"tempfwd_rnaseq.xlsx")
#list2excelout(tempdgeout_rev,temp_lst,"temprev_rnaseq.xlsx")

list2subsettable(tempdgeout,temp_lst,"_temp_all_list.txt")
list2subsettable(tempdgeout_fwd,temp_lst,"_temp_fwd_list.txt")
list2subsettable(tempdgeout_rev,temp_lst,"_temp_rev_list.txt")

#create population files for GO
#each tested set has a different number of testable genes, might affect GO output
list2poptable(nontempdgeout,nontemp_lst,"_nontemp_all_pop.txt")
list2poptable(nontempdgeout_fwd,nontemp_lst,"_nontemp_fwd_pop.txt")
list2poptable(nontempdgeout_rev,nontemp_lst,"_nontemp_rev_pop.txt")
list2poptable(tempdgeout,temp_lst,"_temp_all_pop.txt")
list2poptable(tempdgeout_fwd,temp_lst,"_temp_fwd_pop.txt")
list2poptable(tempdgeout_rev,temp_lst,"_temp_rev_pop.txt")

#compare wt nontemp vs temp controls
wtdgeout = vector("list",(length(1)) )
for (strain in (seq(1:1)) ){
  mutants = grep(temp_lst[[1]],names(pombe_geneOfrac_all),value=T)
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_all),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  wtdgeout[[strain]]=dgplus(pombe_geneOfrac_all,controls,mutants,compare_grp)
  
}

#forward strand
wtdgeout_fwd = vector("list",(length(1)) )
for (strain in (seq(1:1)) ){
  mutants = grep(temp_lst[[1]],names(pombe_geneOfrac_fwd),value=T)
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_fwd),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  wtdgeout_fwd[[strain]]=dgplus(pombe_geneOfrac_fwd,controls,mutants,compare_grp)
  
}

#reverse strand
wtdgeout_rev = vector("list",(length(1)) )
for (strain in (seq(1:1)) ){
  mutants = grep(temp_lst[[1]],names(pombe_geneOfrac_rev),value=T)
  controls = grep(nontemp_lst[[1]],names(pombe_geneOfrac_rev),value=T)
  compare_grp = ifelse(c(controls,mutants) == mutants,1,0 )
  wtdgeout_rev[[strain]]=dgplus(pombe_geneOfrac_rev,controls,mutants,compare_grp)
  
}

#can combine wt output into one excel file for simplicity
wtdgeout_combined = c(wtdgeout,wtdgeout_fwd,wtdgeout_rev)

#output here
wb <- createWorkbook()
datas <- wtdgeout_combined
sheetnames = c("wt_all","wt_forward","wt_reverse")
sheets <- lapply(sheetnames, createSheet, wb = wb)
void <- Map(addDataFrame, datas, sheets)
saveWorkbook(wb, file = "wt_dge_rnaseq.xlsx")


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
