#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-template.tcsh OUTPUT-DIR PARAM-SCRIPT LOOP-BRANCH OBJECT1 OBJECT2
##

if ($#argv != 5) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set object1 = $4
set object2 = $5
set objects = ($object1 $object2)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir winsize"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

# Set parameters
if ($bias_correction == TRUE) then
	set m = "bias"
else
	set m = "nobias"
endif

set main_dir = `echo ${cwd}`
set l1 = "$branch"/"$object1"/loops_unfiltered_"$m"_cpm.tsv.gz
set l2 = "$branch"/"$object2"/loops_unfiltered_"$m"_cpm.tsv.gz
set tss = $genome_dir/tss.bed

# Create uncompressed version of unfiltered loops
echo "Uncompressing $object1 unfiltered loops..." | scripts-send2err
cat $l1 | gunzip >! $outdir/l1.tsv
echo "Uncompressing $object2 unfiltered loops..." | scripts-send2err
cat $l2 | gunzip >! $outdir/l2.tsv

# Soft-filter loops (qcut2 & min_cpm)
echo "Soft-filtering loops (qval <= "$qcut2")..." | scripts-send2err
awk -v q="$qcut2" -v c="$min_cpm" '{if ((NR == 1) || ($7 <= q) && ($5 >= c)){print}}' $outdir/l1.tsv >! $outdir/l1_q2.tsv
awk -v q="$qcut2" -v c="$min_cpm" '{if ((NR == 1) || ($7 <= q) && ($5 >= c)){print}}' $outdir/l2.tsv >! $outdir/l2_q2.tsv

# Generate loops-analysis report
echo "Performing analysis..." | scripts-send2err
Rscript ./code/scripts-loops-diff.r $outdir/l1_q2.tsv $outdir/l2_q2.tsv $object1 $object2 $outdir $winsize $tss $common_log2FC $qcut1 $qcut2 $min_distance $max_distance

# Create bedpe and igv files
echo "Creating bedpe and igv loops files..." | scripts-send2err
mkdir -p $outdir/bedpe_files $outdir/igv_files

#bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/"$object1"_specific_loops.tsv >! $outdir/bedpe_files/"$object1"_specific_loops.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/"$object2"_specific_loops.tsv >! $outdir/bedpe_files/"$object2"_specific_loops.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/common.loops_decreased.tsv >! $outdir/bedpe_files/common.loops_decreased.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/common.loops_increased.tsv >! $outdir/bedpe_files/common.loops_increased.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/common.loops.tsv >! $outdir/bedpe_files/common.loops.bedpe
awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/common.loops_stable.tsv >! $outdir/bedpe_files/common.loops_stable.bedpe

cd $outdir/bedpe_files/ ; wc -l *.bedpe | fgrep -v "total" | awk '{print $1"\t"$2}' > loop_count.tsv
cd $main_dir

#igv
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/common.loops_decreased.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! $outdir/igv_files/common.loops_decreased.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/common.loops_increased.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! $outdir/igv_files/common.loops_increased.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/common.loops.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n > ! $outdir/igv_files/common.loops.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/common.loops_stable.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n > ! $outdir/igv_files/common.loops_stable.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/"$object1"_specific_loops.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! $outdir/igv_files/"$object1"_specific_loops.igv.bed
awk '{if(NR>1) print $1"\t"$2"\t"$4"\t\.\t1\.0"}' $outdir/"$object2"_specific_loops.tsv | sed -e '1itrack graphType=junctions' | sort -k2 -n >! $outdir/igv_files/"$object2"_specific_loops.igv.bed

### APA Analysis ###
if ($APA_diff == TRUE) then
	mkdir -p $outdir/APA
        mkdir -p $outdir/APA/diff
	
	# Perform analysis on the different loops subsets #
	set bedpes =  `cd $outdir/bedpe_files/; ls *.bedpe`

	set job_dir = $outdir/__jdata_APA
	mkdir -p $job_dir

        set hicfile1 = `ls -l inpdirs/tracks/results/tracks.by_*/*/*/"$object1"/filtered.hic | awk '{print $9}' | head -n1`
        set hicfile2 = `ls -l inpdirs/tracks/results/tracks.by_*/*/*/"$object2"/filtered.hic | awk '{print $9}' | head -n1`
  
	set nbed = `ls -l $outdir/bedpe_files/*.bedpe | wc -l`
	
        echo "Computing APA scores on the loop-subsets..." | scripts-send2err
	set jid1 = `sbatch --array=1-$nbed --output="$job_dir/job.%a.out" --error="$job_dir/job.%a.err" ./code/APA.sh $hicfile1 "$bedpes" $APA_resolution $outdir $object1 TRUE diff`
	set jid2 = `sbatch --array=1-$nbed --output="$job_dir/job.%a.out" --error="$job_dir/job.%a.err" ./code/APA.sh $hicfile2 "$bedpes" $APA_resolution $outdir $object2 TRUE diff`

	set jid1 = `echo $jid1 | sed 's/.* //'`
	set jid2 = `echo $jid2 | sed 's/.* //'`

	echo $jid1 >! $job_dir/job1.id
	echo $jid2 >! $job_dir/job2.id

	scripts-send2err "Waiting for job array [$jid1] to complete..."
        scripts-send2err "Waiting for job array [$jid2] to complete..."
	scripts-qsub-wait "$jid1"
       	scripts-qsub-wait "$jid2"
	rm -fr $outdir/APA/diff/*/*/*v*

	echo "Performing in-house APA analysis on the loop-subsets..." | scripts-send2err
	Rscript ./code/scripts-loops-APA-diff.r $outdir/APA/diff/ $APA_resolution $object1 $object2 $URm "$main_dir"/"$outdir"/bedpe_files/loop_count.tsv
endif

# Generate loops sets by quantiles and perform APA analyisis
if ($APA_quantiles == TRUE) then
        mkdir -p $outdir/APA
        mkdir -p $outdir/APA/quantiles
	
	# get quantiles-wise loop subsets	
	if (APA_qfile == ref1) then
		set qf = $outdir/l1_q2.tsv
	else
		awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t0.000000000001\t0.00000000001"}' $outdir/common.loops.tsv > $outdir/common.loops_clean.tsv
            	set qf = $outdir/common.loops_clean.tsv
	endif
	
	awk -v q="$qcut1" -v c="$min_cpm" '{if ((NR == 1) || ($7 <= q) && ($5 >= c)){print}}' $qf >! $outdir/t1.tsv
        awk -v min="$min_distance" -v max="$max_distance" 'NR == 1 || \!/fragment/ && ($4-$2) < max && ($4-$2) > min' $outdir/t1.tsv >! $outdir/t2.tsv        
	awk -v var="$winsize" '{ if ((NR>1)) print $1"\t"($2-var/2)"\t"($2+var/2)"\t"$3"\t"($4-var/2)"\t"($4+var/2)"\t"$5}' $outdir/t2.tsv >! $outdir/l1_q2.bedpe
        rm -f $outdir/t1.tsv $outdir/t2.tsv 

	echo "Obtaining quantile-subsets..." | scripts-send2err
        Rscript ./code/scripts-loops-quantiles.r $outdir/l1_q2.bedpe $outdir/APA/quantiles/

	# perform APA on quantile-subsets
	set job_dir = $outdir/__jdata_APA
	mkdir -p $job_dir

	set hicfile1 = inpdirs/tracks/results/tracks.by_*/*/*/"$object1"/filtered.hic
        set hicfile2 = inpdirs/tracks/results/tracks.by_*/*/*/"$object2"/filtered.hic
        set bedpes = `cd $outdir/APA/quantiles/; ls -l *bedpe | awk '{print $9}'`
	set nbed = `ls -l $outdir/APA/quantiles/*.bedpe | wc -l`

	echo "Computing APA scores on the quantile-subsets..." | scripts-send2err
	set jid1 = `sbatch --array=1-$nbed --output="$job_dir/job.%a.out" --error="$job_dir/job.%a.err" ./code/APA.sh $hicfile1 "$bedpes" $APA_resolution $outdir $object1 TRUE quantiles`
        set jid2 = `sbatch --array=1-$nbed --output="$job_dir/job.%a.out" --error="$job_dir/job.%a.err" ./code/APA.sh $hicfile2 "$bedpes" $APA_resolution $outdir $object2 TRUE quantiles`

        set jid1 = `echo $jid1 | sed 's/.* //'`
        set jid2 = `echo $jid2 | sed 's/.* //'`

        echo $jid1 >! $job_dir/job1.id
        echo $jid2 >! $job_dir/job2.id

        scripts-send2err "Waiting for job array [$jid1] to complete..."
        scripts-send2err "Waiting for job array [$jid2] to complete..."
        scripts-qsub-wait "$jid1"
        scripts-qsub-wait "$jid2"
	rm -fr $outdir/APA/quantiles/*/*/*v* $outdir/APA/quantiles/*bedpe

	# perform in-house analysis
	echo "Performing in-house APA analysis on the quantile-subsets..." | scripts-send2err
	Rscript ./code/scripts-loops-APA-quantiles.r $outdir/APA/quantiles/ $APA_resolution $object1 $object2 $URm

endif

#rm -f $outdir/l*.bedpe $outdir/l*.tsv

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."
