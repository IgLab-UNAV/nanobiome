#Script to read the 16S csv output from EPI2ME, make an OTU table and perform DA analysis.
library(data.table)
library(stringr)
library(edgeR)
library(Glimma)
library(ALDEx2)
# 1. Read 16S files from EPI2ME####
#This includes code from https://github.com/ImagoXV/Epi2me_output_converter
#Files from EPI2ME 16S analysis should be named with the date in format "DD_MM_YYYY.csv"
getwd()
#Read in .csv files. If EPI2ME 16S files are separated by semicolon (;), then use read.csv2:
filenames <- list.files(pattern = "\\.csv")
epi2me_table <- data.frame()
for (i in filenames){
  new_table <- read.csv(i)
  fname <- str_match_all(i, pattern = "(.+).csv")[[1]][2]
  new_table$barcode <- paste(fname, new_table$barcode, sep = ".")
  epi2me_table <- rbind(epi2me_table, new_table)
}
# Remove all rows with "unclassified" barcodes (keep barcodes ending in "barcode__":
retain <- grep(".barcode\\d", epi2me_table$barcode)
epi2me_table.clean <- epi2me_table[retain,]
rm(new_table)
rm(epi2me_table)
# Get a non-redundant vector of all barcodes present in file:
samples <- sort(unique(epi2me_table.clean$barcode))
#Create list object with OTUs in all barcodes:
Multi <- vector(mode='list', length = length(samples))
names(Multi) <- samples
for(i in 1:length(samples)){
  Multi[[i]] <- data.frame(table(epi2me_table.clean[epi2me_table.clean$barcode==samples[i],7]))
  if(ncol(Multi[[i]])==2){
    colnames(Multi[[i]]) <- c("Species", samples[i])
  }
}
#Merge all info for the final OTU counts table: 
mg <- merge(Multi[[1]], Multi[[2]], by = "Species", all = TRUE )
for(i in 3:length(samples)){
  if(ncol(Multi[[i]])==2){
    mg <- merge(mg, Multi[[i]], by = "Species", all = TRUE )
  } else {
    df<- data.frame(mg[,1], NA)
    colnames(df) <- c("Species", samples[i])
    mg <- merge(mg, df, by = "Species", all = TRUE )
  }
}
# Substitute NAs with zero (will also change variables from int to num):
mg[is.na(mg)] <- 0
# See number of counts per column and get rid of barcodes with <2000 counts:
table(colSums(mg[2:ncol(mg)]))
keep.barcodes <- colSums(mg[2:ncol(mg)]) > 2000
mg.clean <- cbind(Species = mg[,1],mg[2:ncol(mg)][keep.barcodes])
colSums(mg.clean[2:ncol(mg.clean)])
#
counts <- mg.clean[,-1]
rownames(counts) <- mg.clean[,1]
#
# Remove rows with < 2 total counts to decrease sparsity:
table(rowSums(counts) > 1)
keep.species <- rowSums(counts) > 1
counts <- counts[keep.species,]
# Remove rows with empty Species name (if any):
table(rownames(counts)=="")
nospecies <- rownames(counts)==""
counts <- counts[!nospecies,]
# 
# Write OTU counts table:
write.table(counts, file = "OTU_table.txt", quote = F, sep = "\t", row.names = T)
#
#
# 2. Differential Abundance (DA) analysis####
# Following recommendations from https://doi.org/10.1038/s41467-022-28034-z, we implement
# two different methods: limma-voom (with TMMwsp normalization) and ALDEx2.
# We do not use rarefied tables (but N.B we have decreased sparsity by removing rows with
# < 2 total counts).
#
#'-------------
# LIMMA-voom
#'-------------
# First we create the metadata dataframe:
sampleinfo <- data.frame(row.names = colnames(counts))
#
# Create column with metadata. CHANGE THIS with your own data (careful here...):
sampleinfo$group <- factor(c(rep("pBIC10",3), rep("pBIC",2),rep("pBI",2), 
                             rep("pBICA",3), "pBIC10", rep("pBIC",2),rep("pBI",2)))
sampleinfo
group <- sampleinfo$group
# Normalize with TMMwsp:
dge <- DGEList(counts, group = group)
dge <- calcNormFactors(dge, method = "TMMwsp")
dge$samples
# Visualizations
barplot(dge$samples$lib.size, names=colnames(dge),las=2, cex.names = 0.6)
log.counts.cpm <- cpm(dge, log = T)
barplot(colSums(log.counts.cpm), names=colnames(log.counts.cpm),las=2, cex.names = 0.6)
boxplot(log.counts.cpm, xlab="", ylab="Log2 counts per million",las=1, cex.axis=0.6, 
        notch = TRUE)
abline(h=median(log.counts.cpm), col="blue")
# Quality check with MDS to check correct clustering, batch effects, etc:
glimmaMDS(dge, width = 920, html = "./MDS.16S.html")
#
# Limma
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
v <- voom(dge, design, plot = TRUE)
head(v$E)
# Now we can check the effect of voom-normalization
par(mfrow=c(1,2))
boxplot(log.counts.cpm, xlab="", ylab="Log2 CPM",las=1, cex.axis = 0.7, main="logCPM")
abline(h=median(log.counts.cpm),col="blue")
boxplot(v$E, xlab="", ylab="Log2 CPM",las=1,cex.axis = 0.7, main="Voom logCPM")
abline(h=median(v$E),col="blue")
# Fit linear model
fit<-lmFit(v)
names(fit)
# Create contrast matrix (for one comparison)
cont.matrix <- makeContrasts(bic_bic10=pBIC - pBIC10, 
                             levels=design)
cont.matrix
# Run statistics between contrats and eBayes shrinkage:
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
# Check number of DA OTUs with different adjusted p-value thresholds
summary(decideTests(fit.cont, p.value = 0.05))
#
# Now get top DA OTUs for the contrast with desired values of adjusted p-value and logFC.
result.limma <- topTable(fit.cont, coef="bic_bic10", p.value = 0.01, lfc = 2, 
                         sort.by="logFC", resort.by = "logFC", number = 20000)
dim(result.limma)
# Check heatmap
daos.limma <- rownames((result.limma))
heatgenes <- (log.counts.cpm[daos.limma,])
pheatmap::pheatmap(heatgenes, scale = "row", fontsize_row = 6, annotation_col = sampleinfo, legend = T)
# write to a file if we wish.
write.table(result.limma, file = "./result.limma.txt", quote = FALSE, sep="\t")
#
#'-----------
# ALDEx2
#'-----------
# Retain only samples from the two groups to be compared and create metadata:
counts.aldex <- counts[,c(1:3,11,4:5,12:13)]
conds <- c(rep("pBIC10", 4), rep("pBIC", 4))
# Run ALDEx2:
x.all <- ALDEx2::aldex(counts.aldex, conds, mc.samples=1000, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
par(mfrow=c(1,2))
ALDEx2::aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
ALDEx2::aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")
result.aldex <- x.all[x.all$we.eBH < 0.1,]
daos.aldex <-rownames(result.aldex)
heatgenes <- (log.counts.cpm[daos.aldex,])
pheatmap::pheatmap(heatgenes, scale = "row", fontsize_row = 6, annotation_col = sampleinfo, legend = T)
#
#Compare results from both tests (most likely all daos.aldex will be included in daos.limma):
intersect(daos.aldex, daos.limma)
#
