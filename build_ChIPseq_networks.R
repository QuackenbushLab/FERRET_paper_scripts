# Next, we build a network that has scores for each TF-gene association.
# We remove files that do not have scores ranging from 1-1000, as defined by UCSC: https://bedtools.readthedocs.io/en/latest/content/general-usage.html
baseRepo <- NULL
fibroblastsConsistentFiles <- c(#"fibroblasts/ATF3/GSE81403_pics_dek007-dek008_peaks_hg19_overlap.bed",
                                #"fibroblasts/MAZ/GSM1003613_hg19_wgEncodeSydhTfbsImr90Mazab85725IggrabPk_hg19_peaks.bed",
                                "fibroblasts/RAD21/GSM935624_hg19_wgEncodeSydhTfbsImr90Rad21IggrabPk_hg19_peaks.bed",
                                #"fibroblasts/RCOR1/GSM1003612_hg19_wgEncodeSydhTfbsImr90Corestsc30189IggrabPk_hg19_peaks.bed",
                                #"fibroblasts/RFX5/GSM1003615_hg19_wgEncodeSydhTfbsImr90Rfx5IggrabPk_hg19_overlap.bed",
                                #"fibroblasts/USF2/GSM2827344_ENCFF212RZM_peaks_GRCh38_hg38_peaks.bed",
                                "fibroblasts/POLR2A/GSM935513_hg19_wgEncodeSydhTfbsImr90Pol2IggrabPk_hg19_peaks.bed")
B_lymphocytesConsistentFiles <- c()#"B_lymphocytes/MAZ/GSM935283_hg19_wgEncodeSydhTfbsGm12878Mazab85725IggmusPk_hg19.bed",
                               #"B_lymphocytes/RFX5/GSM935556_hg19_wgEncodeSydhTfbsGm12878Rfx5200401194IggmusPk_hg19.bed",
                               #"B_lymphocytes/USF2/GSM935558_hg19_wgEncodeSydhTfbsGm12878Usf2IggmusPk_hg19.bed")
fibroblastsLeaveOutFiles <- c("fibroblasts/H2AZ/GSM1003505_hg19_wgEncodeBroadHistoneNhdfadH2azPk_hg19_peaks.bed",
                              "fibroblasts/H2AZ/GSM1003530_hg19_wgEncodeBroadHistoneNhlfH2azPk_hg19_peaks.bed")
                              #"fibroblasts/POLR2A/GSM1330727_MRC5_RNAPII-untreated_DHM025_hg19_peaks.bed",
                              #"fibroblasts/POLR2A/GSM1405132_CS1AN1_N20_vs_input1_peaks_hg19_peaks.bed",
                              #"fibroblasts/POLR2A/GSM1405135_CSB1_N20_vs_input1_peaks_hg19_peaks.bed",
                              #"fibroblasts/POLR2A/GSM1405141_CSB2_N20_vs_input2_peaks_hg19_peaks.bed",
                              #"fibroblasts/POLR2A/GSM1869138_BJ_PolII_MeDiChi-seq_peaks_updatedFormat_hg19.bed",
                              #"fibroblasts/POLR2A/GSM1869150_BJELM_PolII_MeDiChi-seq_peaks_updatedFormat_hg19.bed",
                              #"fibroblasts/POLR2A/GSM2218620_lane7M110p2s2_sequence.merged_hg19_peaks.bed",
                              #"fibroblasts/POLR2A/GSM2218630_C56YUACXX_MULTI_3A_14s008356-1-1_Fousteri_lane5s5p3e8NOUV_sequence.merged_hg19_peaks.bed")
B_lymphocyteLeaveOutFiles <- c(#"B_lymphocytes/ATF3/GSM803508_hg19_wgEncodeHaibTfbsGm12878Atf3Pcr1xPkRep1_hg19_peaks.bed",
                               #"B_lymphocytes/ATF3/GSM803508_hg19_wgEncodeHaibTfbsGm12878Atf3Pcr1xPkRep2_hg19_peaks.bed",
                               "B_lymphocytes/H2AZ/GSM733767_hg19_wgEncodeBroadHistoneGm12878H2azStdPk_hg19.bed",
                               "B_lymphocytes/H2AZ/GSM1003476_hg19_wgEncodeBroadHistoneCd20H2azPk_hg19.bed",
                               "B_lymphocytes/POLR2A/GSM822270_hg19_wgEncodeOpenChromChipGm12878Pol2Pk_hg19.bed",
                               "B_lymphocytes/POLR2A/GSM935608_hg19_wgEncodeSydhTfbsGm12878Pol2s2IggmusPk_hg19.bed",
                               #"B_lymphocytes/POLR2A/GSM1277991_polIIreplicate1_Scaled_BGSub_hg19.bed",
                               #"B_lymphocytes/POLR2A/GSM2295947_un_2_peaks_hg19.bed",
                               "B_lymphocytes/RAD21/GSM803416_hg19_wgEncodeHaibTfbsGm12878Rad21V0416101PkRep1_hg19.bed",
                               "B_lymphocytes/RAD21/GSM803416_hg19_wgEncodeHaibTfbsGm12878Rad21V0416101PkRep2_hg19.bed")
                               #"B_lymphocytes/RAD21/GSM935332_hg19_wgEncodeSydhTfbsGm12878Rad21IggrabPk_hg19.bed")

# Function to combine the files.
CombineFiles <- function(inFiles, tfList, outFile){
  allFilesDF <- do.call(rbind, lapply(1:length(inFiles), function(i){
    # Read BED file.
    file <- inFiles[i]
    bed <- read.table(paste0(baseRepo, file), sep = "\t")

    # Extract the columns of importance.
    tf <- tfList[i]
    gene <- bed[,4]
    score <- bed[,9]
    
    # In some cases, it may be in the 10th column instead.
    if(file %in% c("fibroblasts/POLR2A/GSM1405132_CS1AN1_N20_vs_input1_peaks_hg19_peaks.bed",
                   "fibroblasts/POLR2A/GSM1405135_CSB1_N20_vs_input1_peaks_hg19_peaks.bed",
                   "fibroblasts/POLR2A/GSM1405141_CSB2_N20_vs_input2_peaks_hg19_peaks.bed",
                   "fibroblasts/ATF3/GSE81403_pics_dek007-dek008_peaks_hg19_overlap.bed")){
      score <- bed[,10]
    }else if(file %in% c("fibroblasts/POLR2A/GSM1869138_BJ_PolII_MeDiChi-seq_peaks_updatedFormat_hg19.bed",
                         "fibroblasts/POLR2A/GSM1869150_BJELM_PolII_MeDiChi-seq_peaks_updatedFormat_hg19.bed")){
      score <- bed[,11]
    }else if(file %in% c("fibroblasts/MAZ/GSM1003613_hg19_wgEncodeSydhTfbsImr90Mazab85725IggrabPk_hg19_peaks.bed",
                         "fibroblasts/RAD21/GSM935624_hg19_wgEncodeSydhTfbsImr90Rad21IggrabPk_hg19_peaks.bed",
                         "fibroblasts/RCOR1/GSM1003612_hg19_wgEncodeSydhTfbsImr90Corestsc30189IggrabPk_hg19_peaks.bed",
                         "fibroblasts/RFX5/GSM1003615_hg19_wgEncodeSydhTfbsImr90Rfx5IggrabPk_hg19_overlap.bed",
                         "fibroblasts/USF2/GSM2827344_ENCFF212RZM_peaks_GRCh38_hg38_peaks.bed")){
      score <- bed[,15]
    }
    str(file)
    
    # Remove duplicated genes.
    uniqueGenes <- unique(gene)
    uniqueGeneDF <- do.call(rbind, lapply(uniqueGenes, function(g){
      geneAtIdx <- gene[which(gene == g)]
      scoreAtIdx <- score[which(gene == g)]
      geneWhichMax <- which.max(scoreAtIdx)
      newDF <- data.frame(tf = tf, gene = geneAtIdx[geneWhichMax], score = scoreAtIdx[geneWhichMax])
      return(newDF)
    }))
    
    # Return.
    return(uniqueGeneDF)
  }))
  
  # Consolidate.
  tfGenePairs <- paste(allFilesDF$tf, allFilesDF$gene, sep = "_")
  consolidatedDF <- do.call(rbind, lapply(unique(tfGenePairs), function(pair){
     meanScoreForPair <- mean(allFilesDF[which(tfGenePairs == pair), "score"])
     return(data.frame(tf = allFilesDF[which(tfGenePairs == pair), "tf"][1],
                       gene = allFilesDF[which(tfGenePairs == pair), "gene"][1],
                       score = meanScoreForPair))
  }))
  
  # Save the file.
  str(consolidatedDF)
  write.csv(consolidatedDF, outFile, row.names = FALSE)
}

# Consolidate with leave-one-out.
CombineFiles(inFiles = c(fibroblastsLeaveOutFiles, fibroblastsConsistentFiles),
             tfList = unlist(lapply(c(fibroblastsLeaveOutFiles, fibroblastsConsistentFiles),
                                    function(fileName){
                                      return(strsplit(fileName, "/")[[1]][2])
                                    })),
             outFile = paste0(baseRepo, "/fibroblasts/networkAll.csv"))
CombineFiles(inFiles = c(B_lymphocyteLeaveOutFiles, B_lymphocytesConsistentFiles),
             tfList = unlist(lapply(c(B_lymphocyteLeaveOutFiles, B_lymphocytesConsistentFiles),
                                    function(fileName){
                                      return(strsplit(fileName, "/")[[1]][2])
                                    })),
             outFile = paste0(baseRepo, "/B_lymphocytes/networkAll.csv"))

# Consolidate all permutations.
for(i in 1:length(fibroblastsLeaveOutFiles)){
  str(fibroblastsLeaveOutFiles[-i])
  CombineFiles(inFiles = c(fibroblastsLeaveOutFiles[-i], fibroblastsConsistentFiles),
               tfList = unlist(lapply(c(fibroblastsLeaveOutFiles[-i], fibroblastsConsistentFiles),
                                      function(fileName){
                                        return(strsplit(fileName, "/")[[1]][2])
                                      })),
               outFile = paste0(baseRepo, "/fibroblasts/networkWithout", i, ".csv"))
}
for(i in 1:length(B_lymphocyteLeaveOutFiles)){
  str(B_lymphocyteLeaveOutFiles[-i])
  CombineFiles(inFiles = c(B_lymphocyteLeaveOutFiles[-i], B_lymphocytesConsistentFiles),
               tfList = unlist(lapply(c(B_lymphocyteLeaveOutFiles[-i], fibroblastsConsistentFiles),
                                      function(fileName){
                                        return(strsplit(fileName, "/")[[1]][2])
                                      })),
               outFile = paste0(baseRepo, "/B_lymphocytes/networkWithout", i, ".csv"))
}