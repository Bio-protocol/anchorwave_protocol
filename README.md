# Applying AnchorWave to Address Plant Genome Alignment

## Introduction
anchorwave is a software for  sensitive alignment of genomes with high sequence diversity, extensive structural polymorphism and whole-genome duplication variation



## Installation
Users should first install the following software.

1.	AnchorWave (Song et al., 2022; v1.0.1; https://github.com/baoxingsong/AnchorWave)
2.	samtools (Li et al., 2009; v1.6; http://www.htslib.org)
3.	minimap2 (Li, 2018; v2.17-r941; https://github.com/lh3/minimap2)
4.	ggplot2 (Hadley,2016; v3.3.5; https://ggplot2.tidyverse.org


## Input data
```
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/gff3/zea_mays/Zea_mays.AGPv4.34.gff3.gz
wget http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
```





## Major Steps
### Access the level of chromosome and collineraity
```
samtools faidx Zea_mays.AGPv4.dna.toplevel.fa
samtools faidx Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
less Zea_mays.AGPv4.dna.toplevel.fa.fai 
less Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.fai 

```

### Extract CDS and lift over reference and query genome
```
anchorwave gff2seq -i Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa -o cds.fa
minimap2 -x splice -a -t 10 -p 0.4 -N 20 Zea_mays.AGPv4.dna.toplevel.fa cds.fa > ref.sam
minimap2 -x splice -a -t 10 -p 0.4 -N 20 Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa cds.fa > cds.sam
```


### Visualization the result of collinearity
```
perl alignmentToDotplot.pl Zea_mays.AGPv4.34.gff3 cds.sam > cds.tab
```
and then use R to plot Figure1

### Applying anchorwave to identify the collinearity and WGD
```
anchorwave proali -r Zea_mays.AGPv4.dna.toplevel.fa -i Zea_mays.AGPv4.34.gff3 \
-a cds.sam -as cds.fa -ar ref.sam -s Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -n align1.anchors -R 1 -Q 2 -ns

```
and then use R to plot Figure2

### Base-pair resolved genome alignment
```
anchorwave proali -i Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa -a cds.sam -as cds.fa -ar ref.sam -s \
Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa  -n align.anchors -o align.maf -t 10 -R 1 -Q 2 -f align1.f.maf > log 2>&1
```
The last step outputs alignment result as MAF formation 

