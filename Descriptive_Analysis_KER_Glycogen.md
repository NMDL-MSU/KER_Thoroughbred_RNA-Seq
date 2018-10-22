---
title: RNA-Seq Pipeline Statistics
author: Deborah Velez-Irizarry
date: Mon Oct 22 14:15:20 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Alignment statistics for KER_Glycogen after mapping reads to the reference genome EquCab3  
  
***  
**Code:**  
Parent Directory:  

> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq  
  
Directory/File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;Descriptive_Analysis/Descriptive_Analysis_KER_Glycogen.R  
 
**Input files:**  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;Trimmomatic/trimmomatic_rst.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Condetri/paired/rst_quality_trim_paired.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Condetri/single/rst_quality_trim_single_R1.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Condetri/single/rst_quality_trim_single_R2.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Tophat/summary_alignment.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;Tophat/uniq_depth.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/summary_htseq.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/paired_counts.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/single_counts.txt  
>&nbsp;&nbsp;&nbsp;Cufflinks/MergedGTF/Annotation/annotation.txt
  
**Output files:**  
  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;Descriptive_Analysis_KER_Glycogen/retained_read_stats_KER_Glycogen.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;HTSeq/htseq_counts_KER_Glycogen.txt  
 
***  
### Code  
Clear Environment


```r
rm(list=ls())
```

**Session Information**


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /opt/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas_haswellp-r0.2.20.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitr_1.20
## 
## loaded via a namespace (and not attached):
## [1] compiler_3.5.1  magrittr_1.5    tools_3.5.1     stringi_1.2.3  
## [5] stringr_1.3.1   evaluate_0.10.1
```

**Input Directory**


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/RNA_Seq"
```

### Animal IDs


```r
animID <- paste(rep("G", 40), 1:40, sep="")
animID
```

```
##  [1] "G1"  "G2"  "G3"  "G4"  "G5"  "G6"  "G7"  "G8"  "G9"  "G10" "G11"
## [12] "G12" "G13" "G14" "G15" "G16" "G17" "G18" "G19" "G20" "G21" "G22"
## [23] "G23" "G24" "G25" "G26" "G27" "G28" "G29" "G30" "G31" "G32" "G33"
## [34] "G34" "G35" "G36" "G37" "G38" "G39" "G40"
```

### Load RNA-Seq Pipeline Statistics:  
**Trimmomatic**


```r
adpt <- read.table(paste(dir,"Trimmomatic", "trimmomatic_rst.txt", sep="/"), 
    header=TRUE, row.names=1)
adpt <- adpt[animID,]
dim(adpt)
```

```
## [1] 40  5
```

**Condetri**  
> Paired reads


```r
cond.PE <- read.table(paste(dir,"Condetri/paired", 
    "rst_quality_trim_paired.txt", sep="/"), header=TRUE, row.names=1)
cond.PE <- cond.PE[animID,]
dim(cond.PE)
```

```
## [1] 40  6
```

> Single Reads: R1


```r
cond.SR1 <- read.table(paste(dir,"Condetri/single", 
    "rst_quality_trim_single_R1.txt", sep="/"), header=TRUE, row.names=1)
cond.SR1 <- cond.SR1[paste(animID, "R1", sep="_"),]
dim(cond.SR1)
```

```
## [1] 40  6
```

> Single Reads: R2


```r
cond.SR2 <- read.table(paste(dir,"Condetri/single", 
    "rst_quality_trim_single_R2.txt", sep="/"), header=TRUE, row.names=1)
cond.SR2 <- cond.SR2[paste(animID, "R2", sep="_"),]
dim(cond.SR2)
```

```
## [1] 40  6
```

**Tophat**


```r
tophat <- read.table(paste(dir,"Tophat", "summary_alignment.txt", sep="/"), 
    header=TRUE, row.names=1)
tophat <- tophat[animID,]
dim(tophat)
```

```
## [1] 40  7
```

**Depth**  
> Uniquely mapped reads (includes paired and single reads)


```r
uniq.depth <- read.table(paste(dir,"Tophat", "uniq_depth.txt", sep="/"), 
    header=TRUE, row.names=1)
uniq.depth <- uniq.depth[animID,]
dim(uniq.depth)
```

```
## [1] 40  3
```

> Uniquely mapped reads per chromosome (includes paired and single reads)


```r
uniq.chr.depth <- read.table(paste(dir,"Tophat", "uniq_chr_depth.txt", sep="/"), 
    header=TRUE, row.names=1)
uniq.chr.depth <- uniq.chr.depth[animID,]
dim(uniq.chr.depth)
```

```
## [1] 40 33
```

**HTSeq**  
> Summary HTSeq table


```r
sumHTSeq <- read.table(paste(dir, "HTSeq", "summary_htseq.txt", sep="/"), 
    header=TRUE, row.names=1)
sumHTSeq <- sumHTSeq[animID,]
dim(sumHTSeq)
```

```
## [1] 40  4
```

> Paired read counts table


```r
paired <- read.table(paste(dir, "HTSeq", "paired_counts.txt", sep="/"), 
    header=TRUE, row.names=1)
paired <- paired[,paste(animID, "paired", sep="_")]
dim(paired)
```

```
## [1] 37875    40
```

> Single read count table


```r
single <- read.table(paste(dir, "HTSeq", "single_counts.txt", sep="/"), 
    header=TRUE, row.names=1)
single <- single[,paste(animID, "single", sep="_")]
dim(single)
```

```
## [1] 37875    40
```

> Annotation


```r
annot <- read.table(paste(dir, "Cufflinks/MergedGTF/Annotation", "annotation.txt", sep="/"),
    header=TRUE, row.names=8)
```

 ### Summary Function


```r
summSD <- function(x, dig=4) round(c(summary(x),
     Std.Dev.=sd(x)), dig)[c("Min.", "1st Qu.", "Median", "Mean", 
    "Std.Dev.", "3rd Qu.", "Max.")]
```

### Summary Trimmomatic: Adapter filtering
**Number of Input Read Paires (million reads)**


```r
startPE <- adpt$input.read.pairs*2
names(startPE) <- rownames(adpt)
summSD((startPE/2)/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  27.9128  39.4417  41.5356  43.3955   8.1216  49.1896  59.3701
```

**Number of Input Read (million reads)**


```r
summSD((startPE)/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  55.8257  78.8835  83.0711  86.7910  16.2431  98.3792 118.7403
```

**Number of retained paired reads after filtering adapter sequences (million)**


```r
summSD(adpt$both.surviving*2/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  44.5764  62.3524  66.1126  69.3601  13.2035  79.8916  96.8845
```

Percent retained paired reads


```r
mean((adpt$both.surviving*2) / startPE)*100
```

```
## [1] 79.91442
```

**Number of forward reads surviving without its pair (million)**


```r
summSD(adpt$fwd.only.surviving/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   4.7374   7.5145   8.4837   8.5368   1.7414   9.5337  13.2225
```

**Number of reverse reads surviving without its pair (million)**


```r
summSD(adpt$rev.only.surviving/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   0.0573   0.0797   0.0903   0.0911   0.0147   0.1005   0.1191
```

**Number of single reads surviving (million)**


```r
summSD((adpt$fwd.only.surviving/1e6 + adpt$rev.only.surviving/1e6))
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   4.7946   7.5971   8.5668   8.6278   1.7547   9.6298  13.3340
```

Percent single reads surviving


```r
mean((adpt$fwd.only.surviving + adpt$rev.only.surviving) / startPE)*100
```

```
## [1] 9.941248
```

**Number of retained reads start after adapter trimming**


```r
trimm.out <- rowSums(data.frame((adpt$both.surviving*2), + adpt$fwd.only.surviving + 
    adpt$rev.only.surviving))
names(trimm.out) <- rownames(adpt)
summSD(trimm.out/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  50.2866  70.9942  74.4647  77.9880  14.6825  88.4486 107.6822
```

**Number of dropped reads (million)**


```r
summSD(startPE - trimm.out)/1e6
```

```
##      Min.   1st Qu.    Median      Mean  Std.Dev.   3rd Qu.      Max. 
##  4.897967  7.734613  8.741991  8.803051  1.783959  9.830442 13.573852
```

Percent single reads surviving


```r
mean((startPE - trimm.out) / startPE)*100
```

```
## [1] 10.14433
```

**Percent of retained reads after adapter trimming from total number of sequenced reads**


```r
summSD(trimm.out/startPE)*100
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##    87.44    89.27    89.80    89.86     0.92    90.46    91.54
```

**Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Trimmomatic.  
Output should be zero.


```r
sum(!rowSums(data.frame((adpt$both.surviving*2), 
    (rowSums(data.frame(adpt$fwd.only.surviving, adpt$rev.only.surviving, 
    adpt$dropped)*2)))) == startPE)
```

```
## [1] 0
```

### Summary Condetri: Quality filtering  
**Data Check:** Input reads for Condetri should be equal to output reads for Trimmomatic


```r
start.cond <- rowSums(data.frame(cond.PE$TotReads, cond.SR1$TotReads, 
    cond.SR2$TotReads))
names(start.cond) <- rownames(cond.PE)

# Do the read input for condetri reflect the read output for trimmomatic? 
# Output should be zero.
sum(!start.cond == trimm.out[rownames(start.cond)])
```

```
## [1] 0
```

```r
# Do the row names for paired and single reads condetri output match? 
# Output should be zero.
sum(!rownames(cond.PE) == unlist(lapply(strsplit(rownames(cond.SR1), "_"), 
    function(x) x[[1]][1])) |
    !rownames(cond.PE) == unlist(lapply(strsplit(rownames(cond.SR2), "_"), 
    function(x) x[[1]][1])))
```

```
## [1] 0
```

**Number of paired reads retained after quality trimming (million single)**


```r
summSD((cond.PE$PairReads)/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  28.9048  40.4085  44.6519  46.4033   9.7394  52.1205  67.8049
```

Percent retained paired reads


```r
mean(cond.PE$PairReads / start.cond)*100
```

```
## [1] 59.3862
```

**Number of unpaired reads (million)**


```r
summSD((cond.PE$UnparedReads + cond.SR1$UnparedReads + 
    cond.SR2$UnparedReads)/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  10.7673  15.7566  18.0211  17.8347   3.3555  20.1686  25.1841
```

Percent unpaired reads


```r
mean((cond.PE$UnparedReads + cond.SR1$UnparedReads + 
    cond.SR2$UnparedReads) / start.cond)*100
```

```
## [1] 22.92288
```

**Number of dropped reads after quality filtering (million)**


```r
cond.drop <- data.frame(PE=(cond.PE$TotReads - (cond.PE$PairReads + 
    cond.PE$UnparedReads)),
    SR1=(cond.SR1$TotReads - cond.SR1$UnparedReads),
    SR2=(cond.SR2$TotReads - cond.SR2$UnparedReads))
rownames(cond.drop) <- rownames(cond.PE)
summSD(rowSums(cond.drop/1e6))
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   8.6708  11.5580  13.8492  13.7500   2.5962  15.4183  20.6057
```

Percent dropped


```r
mean(rowSums(cond.drop) / start.cond)*100
```

```
## [1] 17.69092
```

**Number of retained reads after quality filtering (million)**


```r
cond.out <- rowSums(data.frame(cond.PE$PairReads, 
    cond.PE$UnparedReads, cond.SR1$UnparedReads, 
    cond.SR2$UnparedReads))
names(cond.out) <- rownames(cond.PE)
summSD(cond.out/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  40.5519  57.7207  60.9777  64.2380  12.4623  72.3191  89.8720
```

**Percent of retained reads after quality trimming from total number of sequenced reads**


```r
summSD(cond.out/startPE[names(cond.out)])*100
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##    71.11    72.25    73.87    73.97     1.90    75.78    76.82
```

**Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Condetri  
Output should be zero.


```r
sum(!rowSums(data.frame(cond.out, cond.drop)) == start.cond)
```

```
## [1] 0
```

### Summary Tophat: Number of reads aligning to reference genome EquCap3  
**Data Check:** Input reads for Tophat should be equal to output reads for Condreti  


```r
# Do the read input for condetri reflect the read output for trimmomatic? Output should be zero.
start.top <- tophat$total_input_reads
names(start.top) <- rownames(tophat)
sum(sum(!start.top) == cond.out[rownames(tophat)])
```

```
## [1] 0
```

**Number of paired reads aligning to the reference**


```r
summSD(tophat$total_paired_reads/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  24.9221  31.6531  34.3981  36.5423   7.3259  40.1713  51.3904
```

Percent paired reads


```r
mean(tophat$total_paired_reads/start.top)*100
```

```
## [1] 56.94353
```

**Number single reads aligning to the reference**


```r
summSD(tophat$total_unpaired_reads/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   9.4222  12.5320  14.5778  14.2642   2.6225  16.4778  19.5274
```

Paired single


```r
mean(tophat$total_unpaired_reads/start.top)*100
```

```
## [1] 22.36255
```

**Number of unmapped reads**


```r
tophat.drop <- (start.top - (tophat$total_paired_reads + 
    tophat$total_unpaired_reads))
names(tophat.drop) <- rownames(tophat)
summSD(tophat.drop/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   5.2329  10.8001  13.4584  13.4315   3.7128  15.9329  22.2383
```

Percent unmapped


```r
mean(tophat.drop/start.top)*100
```

```
## [1] 20.69392
```

**Total number of mapped reads**


```r
tophat.out <- (tophat$total_paired_reads + tophat$total_unpaired_reads)
names(tophat.out) <- rownames(tophat)
summSD(tophat.out/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  34.9927  45.1848  47.9233  50.8065   9.4286  58.1860  69.4360
```

**Mapping percent from total input reads**


```r
summSD(tophat$read_mapping_rate_percent)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  73.3064  77.3666  78.7738  79.3061   3.2943  81.1514  87.6165
```

Double check you get the same percent of mapped reads


```r
mean(tophat.out/start.top)*100
```

```
## [1] 79.30608
```

**Percent of reads retained (aligned) from total number of sequenced reads**


```r
summSD(tophat.out/startPE[names(tophat.out)])*100
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##    54.27    56.72    58.28    58.65     2.57    60.31    66.32
```

### Samtools: Retain unique reads  
**Number of uniquely alignned reads**


```r
uniq.out <- tophat$total_unique_reads
names(uniq.out) <- rownames(tophat)
summSD(uniq.out/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  34.3730  43.8966  46.8913  49.4232   9.1558  56.5636  67.6416
```

Percent unique reads


```r
mean(uniq.out/tophat.out)*100
```

```
## [1] 97.28592
```

**Number of dropped reads**


```r
summSD(tophat.out - uniq.out)/1e6
```

```
##      Min.   1st Qu.    Median      Mean  Std.Dev.   3rd Qu.      Max. 
## 0.6197310 1.1751725 1.3448190 1.3833034 0.3231297 1.6804188 1.9121410
```

Percent dropped


```r
mean((tophat.out - uniq.out)/tophat.out)*100
```

```
## [1] 2.71408
```

**Number of uniquely algned and properly paired reads**


```r
summSD(tophat$total_unique_properly_paired/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  22.7046  28.3500  30.5686  32.6546   6.4575  36.0813  46.0437
```

**Total number of reads for HTSeq from total number of sequenced reads**


```r
summSD(uniq.out/startPE[names(uniq.out)])*100
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##    52.84    55.11    56.62    57.06     2.61    58.56    65.05
```

**Data Check:** Do the paired, unpaired and dropped reads sum to total input reads for Tophat.  
Output should be zero.


```r
sum(!rowSums(data.frame(tophat.out, tophat.drop)) == start.top)
```

```
## [1] 0
```

```r
### Summary Depth per animal:  
```

Ratio of total nucleotides mapped to total nucleotides sequenced per animal 


```r
depth <- uniq.depth$Average_coverage
names(depth) <- rownames(uniq.depth)
summSD(depth)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  46.1605  55.6121  60.7515  62.7945   9.1862  70.0713  82.5261
```

Summary of total nucleotides mapped per chromosome for EquCap3


```r
ref <- uniq.chr.depth["EquCab3",]
depth.chr <- uniq.chr.depth[-64,]
round(t(apply(depth.chr, 2, summSD)/1e6), 2)
```

```
##         Min. 1st Qu. Median  Mean Std.Dev. 3rd Qu.  Max.
## chr1    2.89    4.05   4.37  4.55     0.91    5.08  6.38
## chr2    1.01    1.27   1.43  1.51     0.32    1.74  2.29
## chr3    0.82    1.00   1.12  1.19     0.25    1.37  1.77
## chr4    0.78    0.97   1.11  1.15     0.22    1.34  1.62
## chr5    1.32    1.57   1.69  1.81     0.36    2.06  2.68
## chr6    2.07    2.62   2.87  3.00     0.61    3.41  4.61
## chr7    1.28    1.53   1.69  1.77     0.34    2.01  2.62
## chr8    0.88    1.12   1.20  1.27     0.23    1.46  1.79
## chr9    0.50    0.66   0.72  0.78     0.17    0.91  1.26
## chr10   1.91    2.64   2.82  3.00     0.61    3.40  4.43
## chr11   2.54    3.28   3.60  3.78     0.78    4.25  5.84
## chr12   1.23    1.90   2.10  2.19     0.49    2.51  3.24
## chr13   1.29    2.52   2.77  2.87     0.71    3.23  4.48
## chr14   1.44    1.83   1.99  2.08     0.40    2.34  2.91
## chr15   0.63    0.80   0.90  0.95     0.18    1.12  1.35
## chr16   0.74    0.95   1.07  1.12     0.22    1.28  1.70
## chr17   0.42    0.54   0.62  0.66     0.15    0.78  1.09
## chr18   1.18    1.47   1.66  1.73     0.34    1.93  2.49
## chr19   0.35    0.44   0.50  0.53     0.11    0.61  0.78
## chr20   0.53    0.67   0.76  0.80     0.18    0.90  1.30
## chr21   0.34    0.44   0.48  0.50     0.09    0.56  0.69
## chr22   0.78    1.16   1.25  1.33     0.30    1.52  1.99
## chr23   0.23    0.29   0.33  0.34     0.07    0.40  0.48
## chr24   0.43    0.54   0.61  0.65     0.14    0.74  1.01
## chr25   0.62    0.77   0.83  0.89     0.18    1.05  1.35
## chr26   0.17    0.22   0.25  0.27     0.06    0.31  0.40
## chr27   1.48    4.09   5.12  5.09     1.53    6.13  9.19
## chr28   0.77    1.04   1.16  1.21     0.26    1.37  1.74
## chr29   0.18    0.22   0.25  0.26     0.05    0.30  0.39
## chr30   0.31    0.42   0.48  0.51     0.14    0.57  1.06
## chr31   0.34    0.57   0.68  0.69     0.19    0.81  1.33
## chrX    0.64    0.80   0.88  0.93     0.19    1.06  1.43
## chrAll 34.37   43.90  46.89 49.42     9.16   56.56 67.64
```

Ratio of total nucleotides mapped to total sequenced per animal


```r
ratio.chr <- do.call(rbind, lapply(colnames(depth.chr), 
    function(x) summSD(depth.chr[,x]/depth.chr[,"chrAll"])* 100))
row.names(ratio.chr) <- colnames(depth.chr)
ratio.chr
```

```
##          Min. 1st Qu. Median   Mean Std.Dev. 3rd Qu.   Max.
## chr1     8.00    8.73   9.25   9.19     0.56    9.56  10.54
## chr2     2.64    2.78   2.94   3.08     0.47    3.10   4.41
## chr3     2.09    2.23   2.35   2.42     0.29    2.49   3.46
## chr4     2.02    2.13   2.25   2.34     0.34    2.39   3.56
## chr5     3.26    3.49   3.60   3.68     0.28    3.78   4.35
## chr6     5.33    5.87   6.04   6.06     0.33    6.28   6.88
## chr7     3.15    3.48   3.59   3.59     0.22    3.71   4.10
## chr8     2.32    2.47   2.55   2.59     0.19    2.66   3.03
## chr9     1.32    1.45   1.55   1.58     0.17    1.65   1.96
## chr10    5.33    5.92   6.04   6.05     0.29    6.22   6.56
## chr11    5.91    7.13   7.50   7.66     0.75    8.21   9.08
## chr12    3.51    4.22   4.41   4.41     0.46    4.66   5.18
## chr13    3.73    5.47   5.95   5.78     0.89    6.39   6.96
## chr14    3.79    4.01   4.18   4.21     0.26    4.34   4.77
## chr15    1.69    1.80   1.87   1.92     0.18    2.01   2.44
## chr16    1.98    2.08   2.20   2.28     0.28    2.37   3.11
## chr17    1.10    1.22   1.28   1.33     0.19    1.36   1.89
## chr18    2.72    3.30   3.53   3.51     0.38    3.68   4.27
## chr19    0.90    0.97   1.04   1.07     0.13    1.12   1.44
## chr20    1.30    1.44   1.57   1.62     0.26    1.68   2.33
## chr21    0.88    0.95   0.99   1.01     0.08    1.05   1.21
## chr22    2.17    2.56   2.71   2.69     0.21    2.83   3.07
## chr23    0.57    0.63   0.67   0.69     0.09    0.71   0.96
## chr24    1.09    1.19   1.27   1.31     0.18    1.36   1.81
## chr25    1.54    1.67   1.75   1.81     0.22    1.87   2.44
## chr26    0.43    0.48   0.51   0.54     0.12    0.54   0.95
## chr27    4.06    9.17  10.49  10.26     2.40   11.67  15.46
## chr28    2.05    2.31   2.39   2.45     0.21    2.58   2.92
## chr29    0.44    0.48   0.51   0.54     0.08    0.55   0.79
## chr30    0.79    0.89   1.01   1.04     0.23    1.11   2.02
## chr31    0.85    1.26   1.39   1.40     0.27    1.47   2.16
## chrX     1.66    1.75   1.81   1.89     0.23    1.91   2.51
## chrAll 100.00  100.00 100.00 100.00     0.00  100.00 100.00
```

### Summary HTSeq:  
**Number of reads processed**


```r
processed.htseq <- rowSums(data.frame(sumHTSeq$Paired_Processed*2, 
    sumHTSeq$Single_Processed))
names(processed.htseq) <- animID
summSD(processed.htseq/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  34.9847  44.8794  47.8438  50.5458   9.3932  57.8057  69.2921
```

**Number of reads processed with no feature**


```r
no.feature.htseq <- rowSums(data.frame(sumHTSeq$Paired_No_Feature*2, 
    sumHTSeq$Single_No_Feature))
names(no.feature.htseq) <- animID
summSD(no.feature.htseq/1e6)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   3.5800   6.4802   7.3737   7.4461   1.7249   8.4345  11.7795
```

**Percent of reads processed with no feature atribute (no counts obtained) per animal**


```r
summSD(no.feature.htseq/processed.htseq*100)
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##   9.6712  13.7955  14.7149  14.6995   1.9701  15.8742  19.7167
```

**Percent of retained processed reads per animal**


```r
summSD(100 - (no.feature.htseq/processed.htseq*100))
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##  80.2833  84.1258  85.2851  85.3005   1.9701  86.2045  90.3288
```

### Summary expressed genes  
 Merge paired and single read counts per animal


```r
anim <- unlist(lapply(strsplit(colnames(paired), "_"), 
    function(x) x[[1]][1]))
counts <- do.call(cbind, lapply(1:ncol(paired), 
    function(x) rowSums(data.frame(paired[,x], single[,x]))))
colnames(counts) <- anim
rownames(counts) <- rownames(paired)
dim(counts)
```

```
## [1] 37875    40
```

**Remove HTSeq statistics from counts (laste five rows)**


```r
tail(counts[,1:5])
```

```
##                             G1      G2      G3      G4      G5
## XLOC_037904               1640     931    1311     294     791
## __no_feature           7305361 4452089 5125900 2802146 3380401
## __ambiguous                  0       0       0       0       0
## __too_low_aQual              0       0       0       0       0
## __not_aligned                0       0       0       0       0
## __alignment_not_unique       0       0       0       0       0
```

```r
counts <- counts[1:(nrow(counts)-5),]
tail(counts[,1:5])
```

```
##               G1  G2   G3  G4  G5
## XLOC_037899    0   0    1   0   0
## XLOC_037900  295 237  349 136 207
## XLOC_037901    1   0    0   3   0
## XLOC_037902   50  96   35 172  31
## XLOC_037903    0   0    1   0   0
## XLOC_037904 1640 931 1311 294 791
```

**Total number of expressed genes (genes with at least one count)**


```r
counts <- counts[rowSums(counts) > 0,]
nrow(counts)
```

```
## [1] 29917
```

**Number of expressed genes per chromosome**


```r
table(as.character(annot[rownames(counts), "chr"]))
```

```
## 
##  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 
##  2317  1554  1572   811   987  1039  1001   970   568   710   614  1609 
## chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chr29  chr3 chr30 
##   960   542   815   549   659   741   374   347   596   338  1267   330 
## chr31  chr4  chr5  chr6  chr7  chr8  chr9  chrX 
##   252  1029  1517  1305  1557  1176   759  1052
```

**Number of autosomal genes**


```r
sum(table(as.character(annot[rownames(counts), "chr"]))[-32])
```

```
## [1] 28865
```

### Summary genes for differential expression analysis   
**Number of gene transcripts with expression counts greater than 2 in all animals**  


```r
# Number of gene transcripts available for downstream analysis:
idx <- apply(counts,1, function(x) sum(x > 2) == ncol(counts))
counts <- counts[idx,]
nrow(counts)
```

```
## [1] 14133
```

Get annotation for expressed genes


```r
annot <- annot[rownames(counts),]
```

**Number of genes per chromosome**


```r
chr.count <- table(as.character(annot$chr))
chr.count
```

```
## 
##  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 
##  1094   730   811   389   548   500   461   503   241   322   290   782 
## chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chr29  chr3 chr30 
##   419   272   377   249   322   367   134   145   275   132   593   140 
## chr31  chr4  chr5  chr6  chr7  chr8  chr9  chrX 
##   113   462   677   608   785   546   349   497
```

**Number of autosomal transcripts**


```r
sum(chr.count[-(length(chr.count))])
```

```
## [1] 13636
```

**Save count HTSeq table**


```r
write.table(counts, paste(dir, "HTSeq", "htseq_counts_KER_Glycogen.txt", sep="/"), 
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
```

### Save RNA-Seq Pipeline read summary


```r
save(startPE, trimm.out, cond.out, tophat.out, uniq.out, depth, 
   depth.chr, ref, ratio.chr, processed.htseq, no.feature.htseq,
   counts, file=paste(getwd(), "retained_read_stats_KER_Glycogen.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
Descriptive_Analysis_KER_Glycogen.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+RNA-Seq Pipeline Statistics
```

