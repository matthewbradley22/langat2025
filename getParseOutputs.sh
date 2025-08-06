#!/bin/bash

for fileNum in {1..8}
do
scp -r mb223@kebnekaise.hpc2n.umu.se:/home/m/mahogny/mystore/dataset/250729_scflavi/fastq/out/\
P36207_100"$fileNum"_S2"$fileNum"_L002_R1_001.out/all-sample/DGE_filtered ~/Documents/ÖverbyLab/FilteredParseOutput/P36207\
_100"$fileNum"_S2"$fileNum"_L002_R1_001

done

#This requires entering password each time through loop
