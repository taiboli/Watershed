library(data.table)
# library(R.utils)

### Example script to make N2 pair index in preparation for Watershed input

all.chrs <- 1:22


for(this.chr in all.chrs){
  MESArv <- fread(paste0("GenotypeMapped/chr", this.chr, "_mapped.txt.gz"))
  names(MESArv) <- c("Ind", "Gene", "Chr", "Start", "End") 
  MESArv$Ind <- unlist(lapply(MESArv$Ind, function(x) strsplit(x, "_")[[1]][[1]]))
  MESArv$Gene <- substr(MESArv$Gene, 1, 15)
  MESArv$RV <- paste0(MESArv$Chr, ":", MESArv$Start)
  MESArv <- MESArv[!duplicated(MESArv),]
  # MESArv$geneRV <- paste0(MESArv$Gene, "+", MESArv$RV)
  setDT(MESArv)
  tmp <- MESArv[, list(RV = paste(RV[order(RV)], collapse="+")), by = list(Gene, Ind)]
  tmp$geneRV <- paste0(tmp$Gene, "+", tmp$RV)
  geneRV.unique <- unique(tmp$geneRV)
  tmp$N2Index <- c(1:length(geneRV.unique))[match(tmp$geneRV, geneRV.unique)]
  # fwrite(tmp, file = paste0("Chr-", this.chr, "-N2PairIndex.txt.gz"),
  #        col.names = T, row.names=F, sep="\t")
  saveRDS(tmp, file = paste0("GenotypeN2PairIndex/Chr", this.chr, "-N2PairIndex.RDS"))
}

