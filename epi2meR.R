#Script to read the 16S csv output from EPI2ME, make an OTU table and analyze with edgeR.
#Includes code from: https://github.com/ImagoXV/Epi2me_output_converter
library(data.table)
library(stringr)
#Files from EPI2ME 16S analysis should be named with the date in format "DD_MM_YYYY.csv"
getwd()
#Read in .csv files. EPI2ME 16S files are separated by semicolon (;), so we use read.csv2:
filenames <- list.files(pattern = "\\.csv")
epi2me_table <- data.frame()
for (i in filenames){
new_table <- read.csv2(i)
fname <- str_match_all(i, pattern = "(.+).csv")[[1]][2]
new_table$barcode <- paste(fname, new_table$barcode, sep = ".")
epi2me_table <- rbind(epi2me_table, new_table)
}
#Remove all rows with unclassified barcodes (keep barcodes ending in "barcode__":
retain <- grep(".barcode\\d", epi2me_table$barcode)
epi2me_table.clean <- epi2me_table[retain,]
rm(new_table)
rm(epi2me_table)
#Get a non-redundant vector of all barcodes present in file:
samples <- sort(unique(epi2me_table.clean$barcode))
#Create list object with species in all barcodes:
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
# See number of counts per column and get rid of barcodes with < 5 counts:
table(colSums(mg[2:ncol(mg)]))
keep.barcodes <- colSums(mg[2:ncol(mg)]) >= 5
mg.clean <- cbind(Species = mg[,1],mg[2:ncol(mg)][keep.barcodes])
# (OPTIONAL) Keep only rows with > 1 count:
table(rowSums(mg[2:ncol(mg)]) > 1)
keep.species <- rowSums(mg[2:ncol(mg)]) > 1
mg.clean <- mg.clean[keep.species,]
# Remove rows with empty Species name (if any):
table(mg.clean[1]=="")
nospecies <- mg.clean[1]==""
mg.clean <- mg.clean[!nospecies,]
# 
# Write OTU counts table:
write.table(mg.clean, file = "OTU_table.txt", quote = F, sep = "\t", row.names = F)
#
#
#DGE-like analysis of differentially represented OTUs with EdgeR####
#
library(edgeR)
library(Glimma)
counts <- mg.clean[,-1]
rownames(counts) <- mg.clean[,1]
#Create metadata dataframe:
sampleinfo <- data.frame(row.names = colnames(counts))
#Create column with metadata (careful here...):
sampleinfo$treatment <- factor(c(rep("cont",4), rep("dis",4),rep("dis",4), rep("cont",4)))
group <- sampleinfo$treatment
dge <- DGEList(counts, group = group)
#Filter rows with less than 5 counts in total
thresh <- counts > 5
table(rowSums(thresh))
keep <- rowSums(thresh) >= 1
summary(keep)
dge.clean <- dge[keep, keep.lib.sizes = FALSE]
# # Alternatively, Filter by Expression with custom args:
# table(filterByExpr(dge, min.count = 3, min.total.count = 8,
#                    large.n = 1, min.prop = 0.25))
# keep <- filterByExpr(dge, min.count = 3, min.total.count = 8,
#                      large.n = 1, min.prop = 0.25)
# dge.clean <- dge[keep, keep.lib.sizes = FALSE]
# 
dge.clean <- calcNormFactors(dge.clean)
dge.clean$samples
# Visualizations
barplot(dge.clean$samples$lib.size, names=colnames(dge),las=2, cex.names = 0.6)
boxplot(dge.clean$counts, xlab="", ylab="Counts",las=1, cex.axis=0.6, 
        notch = TRUE)
logcounts <- cpm(dge.clean, log=T)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=1, cex.axis=0.6, 
        notch = TRUE)
abline(h=median(logcounts), col="blue")
#Heatmap of OTUS with most variation:
var_otus <- apply(logcounts, 1, var)
select_var <- names(sort(var_otus, decreasing=TRUE))[1:100]
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
pheatmap::pheatmap(highly_variable_lcpm, scale = "row", fontsize_row = 3, annotation_col = sampleinfo)
# MDS to check correct clustering, batch effects, etc:
glimmaMDS(dge, width = 920, html = "./MDS.16S.html")
#EdgeR
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design
dge.disp <- estimateDisp(dge.clean, design, robust = TRUE)
sqrt(dge.disp$common.dispersion) #optimal between 0.2 and 0.4
plotBCV(dge.disp)
dge.fit <- glmQLFit(dge.disp, design, robust = TRUE)
dge.fit
colnames(dge.fit$coefficients)
comparison <- makeContrasts(cont-dis, levels = design)
enriched <- glmQLFTest(dge.fit, contrast = comparison)
glimmaVolcano(enriched, dge = dge.clean, main="Cont vs. Dis", html = "./volcano.16S.html")
summary(decideTests(enriched, p.value = 0.1, lfc=1))
enriched.DE <- topTags(enriched, n=500, p.value = 0.1)
dim(enriched.DE)
degs <- abs(enriched.DE$table$logFC) > 1
experiment <- enriched.DE[degs,]
dim(experiment)
# write to a file if we wish.
write.table(experiment, file = "./experiment.txt", quote = FALSE, sep="\t")
 


