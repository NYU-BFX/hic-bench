#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-template.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT(S)
##

if ($#argv != 4) then
  grep '^##' $0
  exit
endif

set outdir = $1
set params = $2
set branch = $3
set objects = ($4)

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$objects" "genome genome_dir bin_size"

# run parameter script
source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------


# echo $inpdirs
# echo $inpdirs[1]
# echo $inpdirs[2]

# 1. hic_matrix
if ( $?inpdirs ) then
  foreach j ( `seq $#inpdirs` )
    set hic_matrix_file_1 = $inpdirs[$j]/filtered.h5
    set hic_matrix_title_1 = $objects[$j]_hicmatrix
    set scaling_factor_for_object = `Rscript ./code/hicseq-collect-library-size.r $library_size_branch $objects[$j]`
    Rscript ./code/hicseq-prepare-ini-hic_matrix.r $hic_matrix_file_1 $hic_matrix_title_1 $scaling_factor_for_object $domains_branch/$objects[$j]/domains.k=001.bed $j $outdir
  end
  # After 1. Combine *hicmatrix.made.ini to one.
  cat $outdir/*hicmatrix.made.ini >> $outdir/combined.01.hic_matrix.ini
endif

# 2. domains-diff
if ( $?domains_diff_branch ) then
  set domains_diff_file_1 = $domains_diff_branch/$object_pair/final_results_forpygenometracks.bed
  set domains_diff_title_1 = Intra_TAD_$object_pair
  Rscript ./code/hicseq-prepare-ini-domains_diff.r $domains_diff_file_1 $domains_diff_title_1 $outdir
endif
# After 2. Combine *domains_diff.made.ini to one.
cat $outdir/*domains_diff.made.ini >> $outdir/combined.02.domains_diff.ini

# 3. compartments
if ( $?compartment_dir ) then
  foreach j ( `seq $#objects` )
    set compartment_file_1 = $compartment_dir/$objects[$j]/pca_activeMarkFix.PC1.bedGraph
    set compartment_title_1 = $objects[$j]_compartments
    Rscript ./code/hicseq-prepare-ini-compartments.r $compartment_file_1 $compartment_title_1 $j $outdir
  end
  # After 3. Combine *compartments.made.ini to one.
  cat $outdir/*compartments.made.ini >> $outdir/combined.03.compartments.ini
endif

# 4. virtual4C

# E1. external_bigwigs
if ( $?external_bigwigs ) then
  foreach j ( `seq $#external_bigwigs` )
    Rscript ./code/hicseq-prepare-ini-external_bigwigs.r $external_bigwigs[$j] `basename $external_bigwigs[$j]` $external_colors[$j] $j $outdir
  end
  # After E1. Combine *-bigwig.made.ini to one.
  cat $outdir/*-bigwig.made.ini >> $outdir/combined.E1.external.data.ini
endif

cat $outdir/combined.*.ini >> $outdir/final.combined.manually.adjust.scales.ini

module load python/cpu/3.6.5

pyGenomeTracks --tracks $outdir/final.combined.manually.adjust.scales.ini --region $inbed --dpi 300 --outFileName $outdir/final.combined.manually.adjust.scales.pdf

# pyGenomeTracks --tracks $outini --region chr15:58,445,797-59,458,364 --outFileName $outdir/region.pdf
# pyGenomeTracks --tracks results/pygenometracks-plot.standard/tracks.by_group.h5.res_5kb/filter.by_sample.mapq_20/align.by_sample.bowtie2/all-samples/external.data.ini --region chr15:58,445,797-59,458,364 --outFileName results/pygenometracks-plot.standard/tracks.by_group.h5.res_5kb/filter.by_sample.mapq_20/align.by_sample.bowtie2/all-samples/region1.pdf
module unload python


if ( [ -z $outdir/final.combined.manually.adjust.scales.pdf ] && [ ! -e $outdir/final.combined.manually.adjust.scales.pdf ] ) then
  rm $outdir/0*.made.ini
endif

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------

# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


