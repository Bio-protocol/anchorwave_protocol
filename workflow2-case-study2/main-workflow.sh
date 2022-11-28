#!/bin/bash

#1.	Download the reference genome and the query genome.
wget https://ftp.ensemblgenomes.org/pub/plants/release-55/gff3/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3.gz
wget https://ftp.ensemblgenomes.org/pub/plants/release-55/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
wget https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz
gunzip *.gz
change format 
sed -i 's/>chr/>/g' Zm-Mo17-REFERENCE-CAU-1.0.fa
#2.	Extract CDS as anchors and lift over to the query and reference genome.
#extract cds 
anchorwave gff2seq -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3 -r Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa  -o cds.fa
#lift over
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zm-Mo17-REFERENCE-CAU-1.0.fa cds.fa > cds.sam
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa cds.fa > ref.sam

# Visualization the result of collinearity and all anchors.
perl alignmentToDotplot.pl Zea_mays.AGPv4.34.gff3 cds.sam > maizecds.tab

#and then use R to plot Figure3 deposited by 2figure3_plot.R

#Applying anchorwave to identify the collinearity
anchorwave genoAli -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3 -as cds.fa -r Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa -a cds.sam -ar ref.sam 
-s Zm-Mo17-REFERENCE-CAU-1.0.fa -n align1.anchor  -IV

#and then use R to plot Figure4 deposited by 3figure4_plot.R

#3.	Perform whole genome alignment.
anchorwave genoAli -i Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.gff3 -as cds.fa -r Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa -a cds.sam -ar ref.sam 
-s Zm-Mo17-REFERENCE-CAU-1.0.fa -n anchors.anchors -o b73tomo17.maf -f  b73tomo17.f.maf -IV
