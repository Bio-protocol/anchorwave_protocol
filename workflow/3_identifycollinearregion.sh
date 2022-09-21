#!/bin/bash

#Applying anchorwave to identify the collinearity and WGD 

anchorwave proali -r Zea_mays.AGPv4.dna.toplevel.fa -i Zea_mays.AGPv4.34.gff3 \
-a cds.sam -as cds.fa -ar ref.sam -s \
Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa -n align1.anchors -R 1 \
-Q 2 -ns
