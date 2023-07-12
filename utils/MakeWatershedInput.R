library(tidyverse)
library(data.table)

### Example script summarizing rare variant calls and genomic annotations to Watershed input format

all.chrs <- 1:22

n2pair.count <- 0 


gencode <- fread("gencode-hg19-TSS-TES.tsv.gz")


annotation.data <- data.table()

for(this.chr in all.chrs){
  
  ## Read genotype data 
  this.chr.rv <- fread(paste0("GenotypeMapped/chr", this.chr, "_mapped.txt.gz"))
  names(this.chr.rv) <- c("IID", "Gene", "Chr", "Start", "End") 
  this.chr.rv$IID <- unlist(lapply(this.chr.rv$IID, function(x) strsplit(x, "_")[[1]][[1]]))
  this.chr.rv$Gene <- substr(this.chr.rv$Gene, 1, 15)
  this.chr.rv$RV <- paste0(this.chr.rv$Chr, ":", this.chr.rv$Start)
  this.chr.rv$geneRV <- paste0(this.chr.rv$Gene, "+", this.chr.rv$RV)
  this.chr.rv$geneInd <- paste0(this.chr.rv$Gene, "+", this.chr.rv$IID)
  
  this.data <- subset(combined.z, geneInd %in% this.chr.rv$geneInd)
  
  ### Match N2 pairs
  tmp <- readRDS(paste0("GenotypeN2PairIndex/Chr", this.chr, "-N2PairIndex.RDS"))
  tmp$Ind <- unlist(lapply(tmp$Ind, function(x) strsplit(x, "_")[[1]][[1]]))
  tmp$geneInd <- paste0(tmp$Gene, "+", tmp$Ind)
  
  this.data$N2 <- tmp$N2Index[match(this.data$geneInd, tmp$geneInd)]
  
  
  ### Set all (gene, individual) pairs with unique RV to have N2 = NA (included in training Watershed)
  count.n2 <- as.data.frame(table(this.data$N2))
  count.n2$Var1 <- as.character(count.n2$Var1)
  
  to.discard <- count.n2$Var1[which(count.n2$Freq < 2)]
  this.data$N2corrected <- ifelse(this.data$N2 %in% to.discard, NA, this.data$N2)
  
  
  ### Discard (gene, individual) pairs where more than 2 pairs have the same RV randomly, and keep rest as N2 pair (evaluation)
  to.discard <- as.integer(count.n2$Var1[which(count.n2$Freq > 2)])
  setDT(this.data)
  setkey(this.data, N2)
  
  to.discard.geneinds <- unlist(lapply(to.discard, function(x){
    this.inds <- this.data[.(x)]
    this.inds <- this.inds$geneInd
    this.redundant <- length(this.inds) - 2
    sample(this.inds, this.redundant)
  }))
  
  
  
  
  
  this.data <- subset(this.data, !geneInd %in% to.discard.geneinds)
  
  
  ### reorder N2 pairs 
  this.data <- this.data[order(this.data$N2corrected),]
  rows.notna <- which(!is.na(this.data$N2corrected))
  this.data$N2 <- NA
  this.data$N2[rows.notna] <- rep(1:(length(rows.notna)/2), each = 2)
  
  this.data$N2 <- this.data$N2 + n2pair.count
  n2pair.count <- n2pair.count + max(this.data$N2, na.rm = T)
  
  this.data[,c("N2corrected") := NULL]
  


  this.chr.rv <- subset(this.chr.rv, geneInd %in% this.data$geneInd)
  this.data <- merge.data.table(this.data[,c("geneInd", "Gene", "IID", 
                                             "RNA", "Splicing", "N2")], 
                                this.chr.rv[,c("geneInd", "RV")], 
                                by = c("geneInd"), all = T)
  this.data$geneRV <- paste0(this.data$Gene, "+", this.data$RV)
  
  
  ### Get VEP annotations 
  this.chr.annotation <- readRDS(paste0("CombinedAnnotations/chr", this.chr, "-combinedAnnotations.RDS"))
  this.data <- merge.data.table(this.data, this.chr.annotation, by = c("geneRV", "RV"))
  
  ### Summarize VEP annotations across genes - here take max for all RVs 
  this.data <- this.data %>%
    group_by(geneInd) %>%
    summarise_all(max)
  
  annotation.data <- rbind(annotation.data, this.data)
  print(this.chr)
}


##############################  Calculation of distTSS, distTES, in-sample MAF, summarized by min

chr.rv.list <- lapply(all.chrs, function(this.chr){
  
  this.chr.rv <- fread(paste0("GenotypeMapped/chr", this.chr, "_mapped.txt.gz"))
  names(this.chr.rv) <- c("IID", "Gene", "Chr", "Start", "End") 
  this.chr.rv$IID <- unlist(lapply(this.chr.rv$IID, function(x) strsplit(x, "_")[[1]][[1]]))
  this.chr.rv$Gene <- substr(this.chr.rv$Gene, 1, 15)
  this.chr.rv$RV <- paste0(this.chr.rv$Chr, ":", this.chr.rv$Start)
  this.chr.rv$geneRV <- paste0(this.chr.rv$Gene, "+", this.chr.rv$RV)
  this.chr.rv$geneInd <- paste0(this.chr.rv$Gene, "+", this.chr.rv$IID)
  
  this.chr.rv <- subset(this.chr.rv, geneInd %in% annotation.data$geneInd)
  
  this.chr.rv.count <- this.chr.rv %>%
    group_by(geneInd) %>%
    summarise(nRV = n())
  setDT(this.chr.rv.count)
  

  this.chr.rv$TSS <- gencode$TSS[match(this.chr.rv$Gene, gencode$Gene)]
  this.chr.rv$TES <- gencode$TES[match(this.chr.rv$Gene, gencode$Gene)]
  this.chr.rv$distTSS <- abs(this.chr.rv$Start - this.chr.rv$TSS)
  this.chr.rv$distTES <- abs(this.chr.rv$Start - this.chr.rv$TES)
  
  this.chr.rv <- this.chr.rv %>%
    group_by(geneInd) %>%
    summarise_all(funs(min), na.rm = T)
  setDT(this.chr.rv)
  this.chr.rv[,c("Start", "End", "TSS", "TES") := NULL]
  this.chr.rv <- merge.data.table(this.chr.rv, this.chr.rv.count, by = c("geneInd"))
  this.chr.rv
})

all.chr.rv.maf.tss.nRV <- rbindlist(chr.rv.list)




### Merge them all together 
### Need to remove RV column

all.chr.rv.maf.tss.nRV[,c("RV", "geneRV"):=NULL]
annotation.data[,c("RV", "geneRV"):=NULL]


common.names <- intersect(names(all.chr.rv.maf.tss.nRV), names(annotation.data))
all.annotations  <- merge.data.table(annotation.data, all.chr.rv.maf.tss.nRV, by = common.names)



all.annotations[,c("geneInd") := NULL]
all.annotations[,c("Chr") := NULL]



### Convert z-scores to p-values
### Here using RNA and Splicing as examples
gene.sign <- sign(all.annotations$RNA)
gene.pvalue <- 2*pnorm(abs(all.annotations$RNA), lower.tail = F)
gene.pvalue <- gene.sign * gene.pvalue
all.annotations$RNA_pvalue <- gene.pvalue

### splicing doesn't have directions 
all.annotations$Splicing_pvalue <- 2*pnorm(abs(all.annotations$Splicing), lower.tail = F)


all.annotations[,c("RNA", "Splicing") := NULL]

### rename columns
names(all.annotations)[1:2] <- c("SubjectID", "GeneName")
all.annotations$N2pair <- all.annotations$N2
all.annotations[,c("N2") := NULL]


fwrite(all.annotations, file = "WatershedInput/AllIndividuals-Training.txt.gz",
       col.names = T, row.names = F, sep = "\t", quote = F)

