# NMDL MSU RNA-Seq Pipeline  
## KER Thoroughbred Project

Brief description of project...

![ker_glycogen](https://user-images.githubusercontent.com/44003875/47320525-ab3c8300-d61f-11e8-9166-53a2cf1ffb80.png)


### Step 1: Filter Adapter Sequences: [Trimmomatic](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/pdf/btu170.pdf/)

**&nbsp;&nbsp;&nbsp;&nbsp;Input:** 

> ID_R1.fastq   
> ID_R2.fastq  

**&nbsp;&nbsp;&nbsp;&nbsp;Output:** 

> ID_R1_PE.fastq  
> ID_R1_SE.fastq  
> ID_R2_PE.fastq  
> ID_R2_SE.fastq 

**Command line to filter adaptor sequences using Trimmomatic**
```bash
> java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE ID_R1.fastq ID_R2.fastq \
> ID_R1_PE.fastq ID_R1_SE.fastq ID_R2_PE.fastq ID_R2_SE.fastq \
> ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

**&nbsp;&nbsp;&nbsp;&nbsp;Parameters:**

> `ILLUMINACLIP`  
> &nbsp;&nbsp;&nbsp;&nbsp;Illumina specific sequences to cut out of sequenced reads (i.e. TruSeq3-PE.fa:2:30:10).
> 
> `LEADING`  and `TRAILING`  
> &nbsp;&nbsp;&nbsp;&nbsp;Will cut basses off the start and end of a read if below a threshold quality. In this example the  
> &nbsp;&nbsp;&nbsp;&nbsp;quality was set to a phred score of 3 to minimize trimming at ends of read. This will be performed  
> &nbsp;&nbsp;&nbsp;&nbsp;latter.
> 
> `SLIDINGWINDOW`  
> &nbsp;&nbsp;&nbsp;&nbsp;Performs a sliding window for trimming and will cut once the average quality falls below a threshold  
> &nbsp;&nbsp;&nbsp;&nbsp;quality. In this example the sliding window is four bases and the quality threshold is set to 15.
>   
> `MINLEN`   
> &nbsp;&nbsp;&nbsp;&nbsp;Minimum length accepted for a read after trimming In this example if a read length if below 36 bases  
> &nbsp;&nbsp;&nbsp;&nbsp;the read is dropped.



### Step 2: Quality Trimming: [Condetri](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0026314&type=printable)

Input  

> 
