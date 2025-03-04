if(!require(rtracklayer)){
  BiocManager::install("rtracklayer")
}
library("rtracklayer")

# For two specific BED files, we fix the format.
GSM1869138_BJ_PolII_MeDiChi-seq_peaks <- NULL
GSM1869138_BJ_PolII_MeDiChi-seq_peaks_out <- NULL
GSM1869150_BJELM_PolII_MeDiChi-seq_peaks <- NULL
GSM1869150_BJELM_PolII_MeDiChi-seq_peaks_out <- NULL
options(scipen = 999999)
FixFormat <- function(inBedFile, outBedFile){
  inBed <- read.table(inBedFile, sep = "\t")
  firstLineRemoved <- inBed[2:nrow(inBed),]
  firstLineRemoved$V2 <- as.numeric(firstLineRemoved$V2)
  firstLineRemoved$V3 <- as.numeric(firstLineRemoved$V3)
  firstLineRemoved$V4 <- as.numeric(firstLineRemoved$V4)
  firstLineRemoved$V5 <- as.numeric(firstLineRemoved$V5)
  firstLineRemoved$V6 <- as.numeric(firstLineRemoved$V6)
  write.table(firstLineRemoved, outBedFile, quote = FALSE, row.names = FALSE, 
              col.names = FALSE, sep = "\t")
}
FixFormat(GSM1869138_BJ_PolII_MeDiChi-seq_peaks, GSM1869138_BJ_PolII_MeDiChi-seq_peaks_out)
FixFormat(GSM1869150_BJELM_PolII_MeDiChi-seq_peaks,GSM1869150_BJELM_PolII_MeDiChi-seq_peaks_out)