#!/bin/bash

#Base-pair resolved genome alignment

anchorwave proali -i Zea_mays.AGPv4.34.gff3 -r Zea_mays.AGPv4.dna.toplevel.fa -a cds.sam -as cds.fa -ar ref.sam -s \
Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa  -n align.anchors -o align.maf -t 10 -R 1 -Q 2 -f align1.f.maf > log 2>&1

