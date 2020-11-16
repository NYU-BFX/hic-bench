#!/bin/bash

input_file=$1
output_file=$2
#input_file='merged_atacPeak147.bed'
#output_file='ESC_atacMerged_tiled_1kb.bed'

awk '($3-$2) < 2000' $input_file | awk '{print $1"\t" sprintf("%.0f",(($3-$2)/2)+$2)}' | awk '{print $1"\t"$2-1"\t"$2+1-1}' > shortCentered.bed
awk '($3-$2) >= 2000' $input_file > long.bed 
bedtools makewindows -b long.bed -w 1 -s 1000 > longTiled.bed
cat shortCentered.bed longTiled.bed | awk '{print $1"\t"$2"\t"$3"\tENH_" ++c"\t1000\t+"}' > $output_file
rm -f shortCentered.bed long.bed longTiled.bed
