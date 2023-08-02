#!/bin/sh

#Download and uncompress genome and GFF file of reference genome
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/fasta/\
zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/gff3/\
zea_mays/Zea_mays.AGPv4.34.gff3.gz
wget http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/\
sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz

#Access the level of chromosome and collineraity
samtools faidx Zea_mays.AGPv4.dna.toplevel.fa
samtools faidx Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa
#Acquire the first contig name "super_16"
less Zea_mays.AGPv4.dna.toplevel.fa.fai 
less Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.fai 


#Extract CDS and lift over reference and query genome
anchorwave gff2seq -i Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa -o cds.fa
#Mapping reference CDS to reference genome
minimap2 -x splice -a -t 10 -k 12 -p 0.4 -N 20 Zea_mays.AGPv4.dna.toplevel.fa cds.fa > ref.sam
#lift over query genome
minimap2 -x splice -a -t 10 -k 12 -p 0.4 -N 20 Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa cds.fa > cds.sam


# Visualization the result of collinearity
perl alignmentToDotplot.pl Zea_mays.AGPv4.34.gff3 cds.sam > cds.tab
