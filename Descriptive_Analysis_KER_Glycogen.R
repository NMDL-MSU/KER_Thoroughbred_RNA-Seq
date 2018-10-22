#' ### Description:  
#' Alignment statistics for KER_Glycogen after mapping reads to the reference genome EquCab3  
#'   
#' ***  
#' **Code:**  
#' Parent Directory:  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
#'   
#' Directory/File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;Descriptive_Analysis/Descriptive_Analysis_KER_Glycogen.R  
#'  
#' **Input files:**  
#' Directory/File:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;Trimmomatic/trimmomatic_rst.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Condetri/paired/rst_quality_trim_paired.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Condetri/single/rst_quality_trim_single_R1.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Condetri/single/rst_quality_trim_single_R2.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Tophat/summary_alignment.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;Tophat/uniq_depth.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/summary_htseq.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/paired_counts.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/single_counts.txt  
#' >&nbsp;&nbsp;&nbsp;Cufflinks/MergedGTF/Annotation/annotation.txt
#'   
#' **Output files:**  
#'   
#' Directory/File:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;Descriptive_Analysis_KER_Glycogen/retained_read_stats_KER_Glycogen.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/htseq_counts_KER_Glycogen.txt  
#'  
#' ***  

#' ### Code  
#' Clear Environment
rm(list=ls())

#' **Session Information**
sessionInfo()

#' **Input Directory**
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"

#' ### Animal IDs
animID <- paste(rep("G", 40), 1:40, sep="")
animID

#' ### Load RNA-Seq Pipeline Statistics:  
#' **Trimmomatic**
adpt <- read.table(paste(dir,"Trimmomatic", "trimmomatic_rst.txt", sep="/"), 
    header=TRUE, row.names=1)
adpt <- adpt[animID,]
dim(adpt)

#' **Condetri**  

#' > Paired reads
cond.PE <- read.table(paste(dir,"Condetri/paired", 
    "rst_quality_trim_paired.txt", sep="/"), header=TRUE, row.names=1)
cond.PE <- cond.PE[animID,]
dim(cond.PE)

#' > Single Reads: R1
cond.SR1 <- read.table(paste(dir,"Condetri/single", 
    "rst_quality_trim_single_R1.txt", sep="/"), header=TRUE, row.names=1)
cond.SR1 <- cond.SR1[paste(animID, "R1", sep="_"),]
dim(cond.SR1)

#' > Single Reads: R2
cond.SR2 <- read.table(paste(dir,"Condetri/single", 
    "rst_quality_trim_single_R2.txt", sep="/"), header=TRUE, row.names=1)
cond.SR2 <- cond.SR2[paste(animID, "R2", sep="_"),]
dim(cond.SR2)

#' **Tophat**
tophat <- read.table(paste(dir,"Tophat", "summary_alignment.txt", sep="/"), 
    header=TRUE, row.names=1)
tophat <- tophat[animID,]
dim(tophat)

#' **Depth**  

#' > Uniquely mapped reads (includes paired and single reads)
uniq.depth <- read.table(paste(dir,"Tophat", "uniq_depth.txt", sep="/"), 
    header=TRUE, row.names=1)
uniq.depth <- uniq.depth[animID,]
dim(uniq.depth)

#' > Uniquely mapped reads per chromosome (includes paired and single reads)
uniq.chr.depth <- read.table(paste(dir,"Tophat", "uniq_chr_depth.txt", sep="/"), 
    header=TRUE, row.names=1)
uniq.chr.depth <- uniq.chr.depth[animID,]
dim(uniq.chr.depth)

#' **HTSeq**  

#' > Summary HTSeq table
sumHTSeq <- read.table(paste(dir, "HTSeq", "summary_htseq.txt", sep="/"), 
    header=TRUE, row.names=1)
sumHTSeq <- sumHTSeq[animID,]
dim(sumHTSeq)

#' > Paired read counts table
paired <- read.table(paste(dir, "HTSeq", "paired_counts.txt", sep="/"), 
    header=TRUE, row.names=1)
paired <- paired[,paste(animID, "paired", sep="_")]
dim(paired)

#' > Single read count table
single <- read.table(paste(dir, "HTSeq", "single_counts.txt", sep="/"), 
    header=TRUE, row.names=1)
single <- single[,paste(animID, "single", sep="_")]
dim(single)

#' > Annotation
annot <- read.table(paste(dir, "Cufflinks/MergedGTF/Annotation", "annotation.txt", sep="/"),
    header=TRUE, row.names=8)

#'  ### Summary Function
summSD <- function(x, dig=4) round(c(summary(x),
     Std.Dev.=sd(x)), dig)[c("Min.", "1st Qu.", "Median", "Mean", 
    "Std.Dev.", "3rd Qu.", "Max.")]


#' ### Summary Trimmomatic: Adapter filtering
#' **Number of Input Read Paires (million reads)**
startPE <- adpt$input.read.pairs*2
names(startPE) <- rownames(adpt)
summSD((startPE/2)/1e6)

#' **Number of Input Read (million reads)**
summSD((startPE)/1e6)

#' **Number of retained paired reads after filtering adapter sequences (million)**
summSD(adpt$both.surviving*2/1e6)

#' Percent retained paired reads
mean((adpt$both.surviving*2) / startPE)*100

#' **Number of forward reads surviving without its pair (million)**
summSD(adpt$fwd.only.surviving/1e6)

#' **Number of reverse reads surviving without its pair (million)**
summSD(adpt$rev.only.surviving/1e6)

#' **Number of single reads surviving (million)**
summSD((adpt$fwd.only.surviving/1e6 + adpt$rev.only.surviving/1e6))

#' Percent single reads surviving
mean((adpt$fwd.only.surviving + adpt$rev.only.surviving) / startPE)*100

#' **Number of retained reads start after adapter trimming**
trimm.out <- rowSums(data.frame((adpt$both.surviving*2), + adpt$fwd.only.surviving + 
    adpt$rev.only.surviving))
names(trimm.out) <- rownames(adpt)
summSD(trimm.out/1e6)

#' **Number of dropped reads (million)**
summSD(startPE - trimm.out)/1e6

#' Percent single reads surviving
mean((startPE - trimm.out) / startPE)*100

#' **Percent of retained reads after adapter trimming from total number of sequenced reads**
summSD(trimm.out/startPE)*100

#' **Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Trimmomatic.  
#' Output should be zero.
sum(!rowSums(data.frame((adpt$both.surviving*2), 
    (rowSums(data.frame(adpt$fwd.only.surviving, adpt$rev.only.surviving, 
    adpt$dropped)*2)))) == startPE)



#' ### Summary Condetri: Quality filtering  
#' **Data Check:** Input reads for Condetri should be equal to output reads for Trimmomatic
start.cond <- rowSums(data.frame(cond.PE$TotReads, cond.SR1$TotReads, 
    cond.SR2$TotReads))
names(start.cond) <- rownames(cond.PE)

# Do the read input for condetri reflect the read output for trimmomatic? 
# Output should be zero.
sum(!start.cond == trimm.out[rownames(start.cond)])

# Do the row names for paired and single reads condetri output match? 
# Output should be zero.
sum(!rownames(cond.PE) == unlist(lapply(strsplit(rownames(cond.SR1), "_"), 
    function(x) x[[1]][1])) |
    !rownames(cond.PE) == unlist(lapply(strsplit(rownames(cond.SR2), "_"), 
    function(x) x[[1]][1])))

#' **Number of paired reads retained after quality trimming (million single)**
summSD((cond.PE$PairReads)/1e6)

#' Percent retained paired reads
mean(cond.PE$PairReads / start.cond)*100

#' **Number of unpaired reads (million)**
summSD((cond.PE$UnparedReads + cond.SR1$UnparedReads + 
    cond.SR2$UnparedReads)/1e6)

#' Percent unpaired reads
mean((cond.PE$UnparedReads + cond.SR1$UnparedReads + 
    cond.SR2$UnparedReads) / start.cond)*100

#' **Number of dropped reads after quality filtering (million)**
cond.drop <- data.frame(PE=(cond.PE$TotReads - (cond.PE$PairReads + 
    cond.PE$UnparedReads)),
    SR1=(cond.SR1$TotReads - cond.SR1$UnparedReads),
    SR2=(cond.SR2$TotReads - cond.SR2$UnparedReads))
rownames(cond.drop) <- rownames(cond.PE)
summSD(rowSums(cond.drop/1e6))

#' Percent dropped
mean(rowSums(cond.drop) / start.cond)*100

#' **Number of retained reads after quality filtering (million)**
cond.out <- rowSums(data.frame(cond.PE$PairReads, 
    cond.PE$UnparedReads, cond.SR1$UnparedReads, 
    cond.SR2$UnparedReads))
names(cond.out) <- rownames(cond.PE)
summSD(cond.out/1e6)

#' **Percent of retained reads after quality trimming from total number of sequenced reads**
summSD(cond.out/startPE[names(cond.out)])*100

#' **Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Condetri  
#' Output should be zero.
sum(!rowSums(data.frame(cond.out, cond.drop)) == start.cond)



#' ### Summary Tophat: Number of reads aligning to reference genome EquCap3  
#' **Data Check:** Input reads for Tophat should be equal to output reads for Condreti  
# Do the read input for condetri reflect the read output for trimmomatic? Output should be zero.
start.top <- tophat$total_input_reads
names(start.top) <- rownames(tophat)
sum(sum(!start.top) == cond.out[rownames(tophat)])

#' **Number of paired reads aligning to the reference**
summSD(tophat$total_paired_reads/1e6)

#' Percent paired reads
mean(tophat$total_paired_reads/start.top)*100

#' **Number single reads aligning to the reference**
summSD(tophat$total_unpaired_reads/1e6)

#' Paired single
mean(tophat$total_unpaired_reads/start.top)*100

#' **Number of unmapped reads**
tophat.drop <- (start.top - (tophat$total_paired_reads + 
    tophat$total_unpaired_reads))
names(tophat.drop) <- rownames(tophat)
summSD(tophat.drop/1e6)

#' Percent unmapped
mean(tophat.drop/start.top)*100

#' **Total number of mapped reads**
tophat.out <- (tophat$total_paired_reads + tophat$total_unpaired_reads)
names(tophat.out) <- rownames(tophat)
summSD(tophat.out/1e6)

#' **Mapping percent from total input reads**
summSD(tophat$read_mapping_rate_percent)

#' Double check you get the same percent of mapped reads
mean(tophat.out/start.top)*100

#' **Percent of reads retained (aligned) from total number of sequenced reads**
summSD(tophat.out/startPE[names(tophat.out)])*100



#' ### Samtools: Retain unique reads  
#' **Number of uniquely alignned reads**
uniq.out <- tophat$total_unique_reads
names(uniq.out) <- rownames(tophat)
summSD(uniq.out/1e6)

#' Percent unique reads
mean(uniq.out/tophat.out)*100

#' **Number of dropped reads**
summSD(tophat.out - uniq.out)/1e6

#' Percent dropped
mean((tophat.out - uniq.out)/tophat.out)*100


#' **Number of uniquely algned and properly paired reads**
summSD(tophat$total_unique_properly_paired/1e6)

#' **Total number of reads for HTSeq from total number of sequenced reads**
summSD(uniq.out/startPE[names(uniq.out)])*100

#' **Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Tophat.  
#' Output should be zero.
sum(!rowSums(data.frame(tophat.out, tophat.drop)) == start.top)



### Summary Depth per animal:  
#' Ratio of total nucleotides mapped to total nucleotides sequenced per animal 
depth <- uniq.depth$Average_coverage
names(depth) <- rownames(uniq.depth)
summSD(depth)

#' Summary of total nucleotides mapped per chromosome for EquCap3
ref <- uniq.chr.depth["EquCab3",]
depth.chr <- uniq.chr.depth[-64,]
round(t(apply(depth.chr, 2, summSD)/1e6), 2)

#' Ratio of total nucleotides mapped to total sequenced per animal
ratio.chr <- do.call(rbind, lapply(colnames(depth.chr), 
    function(x) summSD(depth.chr[,x]/depth.chr[,"chrAll"])* 100))
row.names(ratio.chr) <- colnames(depth.chr)
ratio.chr




#' ### Summary HTSeq:  
#' **Number of reads processed**
processed.htseq <- rowSums(data.frame(sumHTSeq$Paired_Processed*2, 
    sumHTSeq$Single_Processed))
names(processed.htseq) <- animID
summSD(processed.htseq/1e6)

#' **Number of reads processed with no feature**
no.feature.htseq <- rowSums(data.frame(sumHTSeq$Paired_No_Feature*2, 
    sumHTSeq$Single_No_Feature))
names(no.feature.htseq) <- animID
summSD(no.feature.htseq/1e6)

#' **Percent of reads processed with no feature atribute (no counts obtained) per animal**
summSD(no.feature.htseq/processed.htseq*100)

#' **Percent of retained processed reads per animal**
summSD(100 - (no.feature.htseq/processed.htseq*100))


#' ### Summary expressed genes  
#'  Merge paired and single read counts per animal
anim <- unlist(lapply(strsplit(colnames(paired), "_"), 
    function(x) x[[1]][1]))
counts <- do.call(cbind, lapply(1:ncol(paired), 
    function(x) rowSums(data.frame(paired[,x], single[,x]))))
colnames(counts) <- anim
rownames(counts) <- rownames(paired)
dim(counts)

#' **Remove HTSeq statistics from counts (laste five rows)**
tail(counts[,1:5])
counts <- counts[1:(nrow(counts)-5),]
tail(counts[,1:5])


#' **Total number of expressed genes (genes with at least one count)**
counts <- counts[rowSums(counts) > 0,]
nrow(counts)

#' **Number of expressed genes per chromosome**
table(as.character(annot[rownames(counts), "chr"]))

#' **Number of autosomal genes**
sum(table(as.character(annot[rownames(counts), "chr"]))[-32])

#' ### Summary genes for differential expression analysis   
#' **Number of gene transcripts with expression counts greater than 2 in all animals**  
# Number of gene transcripts available for downstream analysis:
idx <- apply(counts,1, function(x) sum(x > 2) == ncol(counts))
counts <- counts[idx,]
nrow(counts)

#' Get annotation for expressed genes
annot <- annot[rownames(counts),]

#' **Number of genes per chromosome**
chr.count <- table(as.character(annot$chr))
chr.count

#' **Number of autosomal transcripts**
sum(chr.count[-(length(chr.count))])


#' **Save count HTSeq table**
write.table(counts, paste(dir, "HTSeq", "htseq_counts_KER_Glycogen.txt", sep="/"), 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")


#' ### Save RNA-Seq Pipeline read summary
save(startPE, trimm.out, cond.out, tophat.out, uniq.out, depth, 
   depth.chr, ref, ratio.chr, processed.htseq, no.feature.htseq,
   counts, file=paste(getwd(), "retained_read_stats_KER_Glycogen.Rdata", sep="/"))



#' ### Run R Script
#+ eval = FALSE
htmlRunR
Descriptive_Analysis_KER_Glycogen.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+RNA-Seq Pipeline Statistics

