#!/bin/bash
main_dir=$1
k27ac=$2
tss=$3

echo 'EC.bg' >> EC_bg.txt
echo 'EX.bg' >> EX_bg.txt
echo 'EC.oe' >> EC_oe.txt
echo 'EX.oe' >> EX_oe.txt

awk 'NR>1' EP_master.tsv > temp.txt

awk '{				
if ($14 > 0)			
	print 1/$14		
else				
	print $14		
}' temp.txt >> EC_bg.txt


awk '{				
if ($8 > 0)			
	print 1/$8		
else				
	print $8		
}' temp.txt >> EX_bg.txt


paste EP_master.tsv EC_bg.txt EX_bg.txt > temp2.txt
awk 'NR>1' temp2.txt > temp3.txt


awk '{				
if ($21 > 0)			
	print $4/$21		
else				
	print "NA"		
}' temp3.txt >> EC_oe.txt


awk '{				
if ($22 > 0)			
	print $5/$22		
else				
	print "NA"		
}' temp3.txt >> EX_oe.txt


sort -k4,4b $main_dir/$tss | cut -f 2-6 > tss_sort.bed
sort -k5,5b $main_dir/$k27ac | cut -f 2-5 > k27ac_sort.bed
paste temp2.txt EC_oe.txt EX_oe.txt | sort -k2,2b > temp4.txt
join -1 2 -2 3 temp4.txt tss_sort.bed | sort -k3,3b > temp5.txt
join -1 3 -2 4 temp5.txt k27ac_sort.bed > temp6.txt

echo "EP.id id1 id2 EC EX EC.EX_sum EC.EX_diff e.hub e.interactivity e.impact p.EX.sum e.phub e.p.interactivity p.hub pp.hub pe.hub pp.interactivity pe.interactivity p.interactivity chr start end distance loop.activity EC.bg EX.bg EC.oe	 EX.oe tss_start tss_end ensembleID strand enh_start enh_end enh_signal"  >> EP_master_plus.tsv   
awk -v OFS="\t" '{print $3,$2,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35}' temp6.txt >> EP_master_plus.tsv

rm -f EC_bg.txt EX_bg.txt EC_oe.txt EX_oe.txt temp*.txt tss_sort.bed k27ac_sort.bed
