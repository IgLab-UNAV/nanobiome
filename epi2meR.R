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
