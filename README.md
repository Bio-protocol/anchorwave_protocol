# Applying AnchorWave to Address Plant Genome Alignment

## Introduction
anchorwave is a software for  sensitive alignment of genomes with high sequence diversity, extensive structural polymorphism and whole-genome duplication variation.

The process generally includes three steps: 1.	 Extract CDS as anchors, 2. lift over to the query and reference genome, 3. Perform whole genome alignment and the other steps that visualizes the relationship of two genomes.

There are two functions implemented in AnchorWave for genome alignment, “genoAli” and “proali”. “genoAli” is suitable for genome alignment without translocation or chromosome fusion. “genoAli” is design to alignment the genomes from different accessions of the same species. It’s also recommended that “genoAli” is used at related species with few structural variations. “proali” is suitable for genome alignment with translocation variation, chromosome fusion or even whole genome duplication.

In this case study, we use ["proali"](https://github.com/Bio-protocol/anchorwave_protocol/tree/master/workflow) for Genome alignment with relocation variation, chromosome fusion or whole genome duplication.
In the case study 2, we use ["genoAli"](https://github.com/Bio-protocol/anchorwave_protocol/tree/master/workflow2-case-study2) for Genome alignment without translocation rearrangement while with inversions.
## Installation
Users should first install the following software.

1.	AnchorWave (Song et al., 2022; v1.0.1; https://github.com/baoxingsong/AnchorWave)
2.	samtools (Li et al., 2009; v1.6; http://www.htslib.org)
3.	minimap2 (Li, 2018; v2.17-r941; https://github.com/lh3/minimap2)
4.	ggplot2 (Hadley,2016; v3.3.5; https://ggplot2.tidyverse.org

## Installation AnchorWave
```
git clone https://github.com/baoxingsong/anchorwave.git
cd anchorwave
cmake ./
make
```
```
#Installation using conda
conda install -c bioconda -c conda-forge anchorwave
```




## Major Steps
## Case study 1:Align two genomes with relocation variation, chromosome fusion or whole genome duplication variation
In the first case study, we align the sorghum genome to the maize genome via AnchorWave for a step-by-step illustration.
### 1.Input data
```
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/gff3/zea_mays/Zea_mays.AGPv4.34.gff3.gz
wget http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
gunzip *.gz
```

#### Access the level of chromosome and collineraity

if we are familar with the genomes, the step can skip.

```
#Access the level of chromosome and collineraity
samtools faidx Zea_mays.AGPv4.dna.toplevel.fa
samtools faidx Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
#Acquire the first contig name "super_16" 
less Zea_mays.AGPv4.dna.toplevel.fa.fai 
less Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.fai 

```

### 2.Extract CDS and lift over reference and query genome

For AnchorWave, “gff2seq” function is used for extracting CDS sequences, ‘r’, ‘i’ and ‘o’ represents input reference FASTA, GFF(3) file and output files. For minimap2, “x splice” function represents long-read splice alignment, ‘t’ represents the number of threads, ‘k’ represents k-mer size, ‘p’ represents min score ratio and ‘a’ represents that output file is SAM format. Other parameters are default. 

```
## extract cds from maize ref ##
anchorwave gff2seq -i Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa -o cds.fa
## lift over ref ##
minimap2 -x splice -a -t 10 -p 0.4 -N 20 Zea_mays.AGPv4.dna.toplevel.fa cds.fa > ref.sam
## lift over query ##
minimap2 -x splice -a -t 10 -p 0.4 -N 20 Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa cds.fa > cds.sam
```


#### Visualization the result of collinearity
At this step, we visually check the whole genome duplication difference and chromosome rearrangements between the reference genome and query genome.
[alignmentToDotplot.pl](https://github.com/Bio-protocol/anchorwave_protocol/blob/master/workflow/alignmentToDotplot.pl) transforms sam file into a plain file
```
#convert to the SAM file
perl alignmentToDotplot.pl Zea_mays.AGPv4.34.gff3 cds.sam > cds.tab
```
and then use R to plot Figure1.
Please noted that we need to modify the chromosomes number manually if we change the species.
```
#Use this piece of R to draw a dotplot.
library(ggplot2)
#Transform Coordinates using a function.
changetoM <- function ( position ){
position=position/1000000;
paste(position, "M", sep="")}
#Read gene position, belong to which chromosome and so on
data =read.table("cds.tab")
#Select all euchromosomes as factor.
#modify the chromosomes number manually if we change the species.
data = data[which(data$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data = data[which(data$V3 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data$V1 = factor(data$V1, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
data$V3 = factor(data$V3, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
#Using ggplot2 to plot a dotplot and beautify it.
#modify the chromosomes number manually if we change the species.
png("E:\\R_Scripts\\figure1.png")
ggplot(data=data, aes(x=V4, y=V2)) +geom_point(size=0.5, aes(color=V5)) +
facet_grid(V1 ~ V3, scales="free",space="free") +labs(x="sorghum", y="maize")+
scale_x_continuous(labels=changetoM) + scale_y_continuous(labels=changetoM) +
theme(axis.line = element_blank(),
panel.background = element_blank(),
panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
axis.text.y = element_text( colour = "black"),
legend.position='none',
axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
dev.off()

```
![image](https://github.com/Bio-protocol/anchorwave_protocol/blob/master/graphs/figure1.png)

Figure 1. The identified all anchors between the maize and sorghum genome.

#### Applying anchorwave to identify the collinearity and WGD

For AnchorWave, "proali" function is used for whole genome alignment, ‘as’ represents anchor sequence files, ‘i’and ‘r’ represents input reference GFF(3) and FASTA file . ‘a’ and “ar” represents SAM file generated by mapping conserved sequence to reference and query genome, ‘s’ represents query genome file, ‘n’ represents output anchor file.“-R 1 -Q 2” represents the number of WGD.

```
anchorwave proali -r Zea_mays.AGPv4.dna.toplevel.fa -i Zea_mays.AGPv4.34.gff3 \
-a cds.sam -as cds.fa -ar ref.sam -s Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -n align1.anchors -R 1 -Q 2 -ns

```
and then use R to plot Figure2

```
changetoM <- function ( position ){
position=position/1000000;
paste(position, "M", sep="")}
library(ggplot2)
data =read.table("align1.anchors", header=TRUE)
#Select all euchromosomes as factor
#modify the chromosomes number manually if we change the species.
data = data[which(data$refChr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data = data[which(data$queryChr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data$refChr = factor(data$refChr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
data$queryCh = factor(data$queryChr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
png("E:\\R_Scripts\\figure2.png")
ggplot(data=data, aes(x=queryStart, y=referenceStart))+
geom_point(size=0.5, aes(color=strand)) +
facet_grid(refChr~queryChr, scales="free", space="free") +
labs(x="sorghum", y="maize")+scale_x_continuous(labels=changetoM) +
scale_y_continuous(labels=changetoM) +
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill =NA,color="black", size=0.5, linetype="solid"),
        axis.text.y = element_text( colour = "black"),
        legend.position='none',
        axis.text.x = element_text(angle=300,hjust=0, vjust=0.5, colour = "black") )
dev.off()
```
![image](https://github.com/Bio-protocol/anchorwave_protocol/blob/master/graphs/figure2.png)

Figure 2. The identified collinear anchors between the maize and sorghum genome.
### 3.Base-pair resolved genome alignment
For AnchorWave, "proali" function is used for whole genome alignment, ‘as’ represents anchor sequence files, ‘i’, ‘r’ and ‘o’ represents input reference GFF(3), FASTA file and output MAF file. ‘a’ and “ar” represents SAM file generated by mapping conserved sequence to reference and query genome, ‘s’ represents query genome file, ‘n’ represents output anchor file.“-R 1 -Q 2” represents the number of WGD.  We perform this step using a workstation with random access memory (RAM) of 128 Gb at least. and the step spent about 36 hours and more.

```
anchorwave proali -i Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa -a cds.sam -as cds.fa -ar ref.sam -s \
Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa  -n align.anchors -o align.maf -t 1 -R 1 -Q 2 -f align1.f.maf > log 2>&1
```

## Case study 2: Align two maize genomes without translocation rearrangement while with inversions

The process uses function "genoAli" that is suitable for genome alignment without translocation or chromosome fusion.we align different maize genomes between B73 and Mo17 to illustrate “genoAli” function. 

### 1.Input data
```
wget https://ftp.ensemblgenomes.org/pub/plants/release-55/gff3/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3.gz
wget https://ftp.ensemblgenomes.org/pub/plants/release-55/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
wget https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz
gunzip *.gz
#change format 
sed -i 's/>chr/>/g' Zm-Mo17-REFERENCE-CAU-1.0.fa
```



### 2.Extract CDS and lift over reference and query genome

For AnchorWave, “gff2seq” function is used for extracting CDS sequences, ‘r’, ‘i’ and ‘o’ represents input reference FASTA, GFF(3) file and output files. For minimap2, “x splice” function represents long-read splice alignment, ‘t’ represents the number of threads, ‘k’ represents k-mer size, ‘p’ represents min score ratio and ‘a’ represents that output file is SAM format. Other parameters are default. 

```
#extract cds from maize ref 
anchorwave gff2seq -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3 -r Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa  -o cds.fa
# lift over re
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zm-Mo17-REFERENCE-CAU-1.0.fa cds.fa > cds.sam
# lift over quer
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa cds.fa > ref.sam

```
We can use R to visualize the full-length CDS mapping result.

#### Visualization the result of collinearity

```
#convert to the SAM file
perl alignmentToDotplot.pl Zea_mays.AGPv4.34.gff3 cds.sam > maizecds.tab
```
and then use R to plot Figure1
```
#Use this piece of R to draw a dotplot.
library(ggplot2)
#Transform Coordinates using a function.
changetoM <- function ( position ){
  position=position/1000000;
  paste(position, "M", sep="")
}
#Read gene position, belong to which chromosome and so on
data =read.table("maizecds.tab")
#Select all euchromosomes as factor, modify the chromosomes number manually.
data = data[which(data$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10" )),]
data = data[which(data$V3 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data$V1 = factor(data$V1, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10" ))
data$V3 = factor(data$V3, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10" ))
#Using ggplot2 to plot a dotplot and beautify it.
png("E:\\plot\\gaoliang\\figure3.png")
ggplot(data=data, aes(x=V4, y=V2))+geom_point(size=0.5, aes(color=V5))+facet_grid(V1~V3, scales="free", space="free" ) +
  labs(x="Mo17", y="B73")+scale_x_continuous(labels=changetoM) + scale_y_continuous(labels=changetoM) +
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.y = element_text( colour = "black"),
        legend.position='none',
        axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black") )
dev.off()

```
![image](https://github.com/Bio-protocol/anchorwave_protocol/blob/master/graphs/figure3.png)

Figure 3. The identified all anchors between the maize B73 and maize Mo17 genome.

#### Applying anchorwave to identify the collinearity and WGD
In this step, We use the “genoAli” function to identify collinear region and use R package ggplot2 to plot the result. Draw the graph with start location of query anchors as x-axis and start location of reference anchors as y-axis.
For	AnchorWave, “genoAli” function is used for whole genome alignment without WGD, ‘as’ represents anchor sequence files, ‘i’, ‘r’ and ‘o’ represents input reference GFF(3), FASTA file and output MAF file. ‘a’ and “ar” represents SAM file generated by mapping conserved sequence to reference and query genome, ‘s’ represents query genome file, ‘n’ represents output anchor files. “IV” represents calling inversion. Other parameters are default.

```
anchorwave genoAli -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3 -as cds.fa -r Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa -a cds.sam -ar ref.sam -s Zm-Mo17-REFERENCE-CAU-1.0.fa -n align1.anchors -IV

```
and then use R to plot Figure2

```
changetoM <- function ( position ){
position=position/1000000;
paste(position, "M", sep="")}
library(ggplot2)
data =read.table("align1.anchors", header=TRUE)
#Select all euchromosomes as factor
data = data[which(data$refChr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data = data[which(data$queryChr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data$refChr = factor(data$refChr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
data$queryCh = factor(data$queryChr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
png("E:\\R_Scripts\\figure2.png")
ggplot(data=data, aes(x=queryStart, y=referenceStart))+
geom_point(size=0.5, aes(color=strand)) +
facet_grid(refChr~queryChr, scales="free", space="free") +
labs(x="sorghum", y="maize")+scale_x_continuous(labels=changetoM) +
scale_y_continuous(labels=changetoM) +
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill =NA,color="black", size=0.5, linetype="solid"),
        axis.text.y = element_text( colour = "black"),
        legend.position='none',
        axis.text.x = element_text(angle=300,hjust=0, vjust=0.5, colour = "black") )
dev.off()
```
![image](https://github.com/Bio-protocol/anchorwave_protocol/blob/master/graphs/figure4.png)

Figure 4. The identified collinear anchors between the maize B73 and maize Mo17 genome.
### 3.Base-pair resolved genome alignment
For	AnchorWave, “genoAli” function is used for whole genome alignment without WGD, ‘as’ represents anchor sequence files, ‘i’, ‘r’ and ‘o’ represents input reference GFF(3), FASTA file and output MAF file. ‘a’ and “ar” represents SAM file generated by mapping conserved sequence to reference and query genome, ‘s’ represents query genome file, ‘n’ represents output anchor files. “IV” represents calling inversion. Other parameters are default. We perform this step using a workstation with random access memory (RAM) of 128 Gb at least and the step spent about 36 hours and more


```
anchorwave genoAli -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3 -as cds.fa -r Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa -a cds.sam -ar ref.sam 
-s Zm-Mo17-REFERENCE-CAU-1.0.fa -n anchors.anchors -o b73tomo17.maf -f  b73tomo17.f.maf -IV
```

## Expecteed results

For case study 1:The last step outputs alignment result as [MAF formation](https://github.com/Bio-protocol/anchorwave_protocol/blob/master/output/align_part.txt) and the out alignment file is about 5.2GB.
For case study 2:The last step also outputs alignment result as MAF formation and the out alignment file is about 5.5GB.




