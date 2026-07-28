if(!require(rtracklayer)){
  BiocManager::install("rtracklayer")
}
library("rtracklayer")
hg19Loc <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Cistrome_Data_Browser/genomes/hg19.gtf"
hg38Loc <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Cistrome_Data_Browser/genomes/hg38.gtf"
hg19BedLoc <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Cistrome_Data_Browser/genomes/hg19.bed"
hg38BedLoc <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/Cistrome_Data_Browser/genomes/hg38.bed"
makeTSSBed <- function(fileLoc, bedFileLoc){
  # Read the GTF file.
  gtf <- rtracklayer::import(fileLoc)
  gtf_df <- as.data.frame(gtf)

  # Filter to include only the protein coding "gene" annotations. Exclude
  # the Y chromosome and mitochondrial DNA.
  # Find out what readthrough_gene and overlapping_locus mean.
  gtf_df$seqnames <- as.character(gtf_df$seqnames)
  gtf_df$strand <- as.character(gtf_df$strand)
  gtf_df$gene_name <- as.character(gtf_df$gene_name)
  gtf_df$type <- as.character(gtf_df$type)
  proteinCoding <- gtf_df[which(gtf_df$gene_type == "protein_coding"),]
  noY <- proteinCoding[which(proteinCoding$seqnames != "chrY"),]
  noM <- noY[which(noY$seqnames != "chrM"),]
  onlyGenes <- noM[which(noM$type == "gene"),]
  
  # Obtain the promoter region for each gene.
  promoterList <- lapply(1:nrow(onlyGenes), function(i){
    
    # If it is a negative strand, promoter is -250, +750 bp.
    # If it is a positive strand, promoter is -750, +250 bp.
    promoterStart <- 0
    promoterEnd <- 0
    if(onlyGenes[i, "strand"] == "+"){
      promoterStart <- onlyGenes[i, "start"] - 750
      promoterEnd <- onlyGenes[i, "start"] + 250
    }else{
      promoterEnd <- onlyGenes[i, "end"] + 750
      promoterStart <- onlyGenes[i, "end"] - 250
    }
    
    # Add chromosome, start, end, and gene name.
    promoterDF <- data.frame(chromosome = onlyGenes[i, "seqnames"], start = promoterStart,
                             end = promoterEnd, name = onlyGenes[i, "gene_name"])

    # Return a data frame.
    return(promoterDF)
  })
  
  # Concatenate all data frames to obtain the final BED file.
  promoters <- do.call(rbind, promoterList)
  str(promoters)
  
  # Save the results in BED format.
  write.table(promoters, bedFileLoc, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
}
makeTSSBed(hg19Loc, hg19BedLoc)
makeTSSBed(hg38Loc, hg38BedLoc)