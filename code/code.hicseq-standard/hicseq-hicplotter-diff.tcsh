#!/bin/tcsh
source ./code/code.main/custom-tcshrc     # shell settings

##
## USAGE: hicseq-hicplotter-diff.tcsh OUTPUT-DIR PARAM-SCRIPT BRANCH OBJECT1 OBJECT2
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

# read variables from input branch
source ./code/code.main/scripts-read-job-vars $branch "$object1 $object2" "genome genome_dir bin_size"

# run parameter script
#source $params

# create path
scripts-create-path $outdir/

# -------------------------------------
# -----  MAIN CODE BELOW --------------
# -------------------------------------

if ($object1 == $object2) goto goto_done

# Create the working directory
if ($?TMP) then
  set tempdir = $TMP
else
  set tempdir = $outdir
endif
set workdir = $tempdir/work
mkdir -p $workdir

# run parameter script
source $params

# create input matrices for HiCplotter
set chrom_included = `echo $regions | tr ' ' '\n' | cut -d':' -f1 | sort -u | tr '\n' '|' | sed 's/|$//'`
foreach f (`cd $branch/$object1; ls -1 matrix.*.tsv matrix.*.RData | grep -vwE "$chrom_excluded" | grep -wE "$chrom_included"`)
  set chr = `echo $f | cut -d'.' -f2`
  scripts-send2err "Processing matrix $f..."
  if ((-e $branch/$object1/$f) && (-e $branch/$object2/$f)) then
    # extract matrices from RData file
    if (`echo $branch/$object1/$f | grep -c '\.RData$'` == 1 & `echo $branch/$object2/$f | grep -c '\.RData$'`) then
      Rscript ./code/hic-matrix.r matrices -v -o $workdir/tmp $branch/$object1/$f
      Rscript ./code/hic-matrix.r matrices -v -o $workdir/tmp $branch/$object2/$f
    else
      mkdir -p $workdir/tmp
      cat $branch/$object1/$f >! $workdir/tmp/matrix.$object1.k=001.tsv
      cat $branch/$object2/$f >! $workdir/tmp/matrix.$object2.k=001.tsv
    endif

    # convert matrices and run hicplotter
    scripts-send2err "Converting matrices into hicplotter format..."
    foreach mat ($workdir/tmp/matrix.*.tsv)
      set pref = `basename $mat .tsv | sed 's/^matrix\.//'`.$chr
      set inpmat = $pref.matrix.txt
      Rscript ./code/create-hicplotter-matrix.r $workdir/$inpmat $mat 1
    end
    rm -rf $workdir/tmp
  endif
end


# Run hicplotter, each region separately
foreach region ($regions)
  echo $region
  scripts-send2err "region = $region"
  set chrom = `echo $region | cut -d':' -f1`
  set start = `echo $region | cut -d'-' -f1 | cut -d':' -f2`
  set stop = `echo $region | cut -d'-' -f2`
  set start_bin = `echo $start/$bin_size | bc`
  set stop_bin = `echo $stop/$bin_size | bc`
  set hic_matrices = ($object2.k=001.$chrom.matrix.txt $object1.k=001.$chrom.matrix.txt)
  set hic_matrix_no = `echo $hic_matrices | tr ' ' '\n' | wc -l`

  echo $chrom
  echo $start $stop
  echo $start_bin $stop_bin
  echo $hic_matrices
  echo $hic_matrix_no
  #set tiles2 = ()
  #foreach t ($tiles)
  #  set tiles2 = ($tiles2 `readlink -f $t`)
  #end
  #set tiles_csv = `echo $tiles2 | sed 's/ /,/g'`
  #set tiles_labels_csv = `echo $tiles_labels | sed 's/ /,/g'`

  #Get the region labels
  #set region_labels = ()
  #Get the loop beds
  #set loop_beds = ()
  #Get all bedgraphs
  set bedgraph_files = ()
  set domain_files = `echo $beds`
  set sample_labels = "$object2 $object1"
  #Get all bedgraph labels
  set bedgraph_files_labels = ()

  foreach i ( `seq 1 1 $hic_matrix_no` )
    #set region_labels = ( $region_labels $region )
    #set loop_beds = ( $loop_beds $loop_bed )
    set bedgraph_files = ( $bedgraph_files $bedgraphs )
    set bedgraph_files_labels = ($bedgraph_files_labels $bedgraph_labels)
  end

  #Run HiC-Plotter
  set p = `pwd`
  if ($highlight == 1) then
    set highlight_opt = "-high $highlight -hf `readlink -f $highlight_bed`"
  else 
    set highlight_opt = 0 
  endif
  set hicplotter_abs_path = `readlink -f $hicplotter_path`
  cd $workdir
  foreach filetype (pdf png)
    echo python $hicplotter_abs_path -v -f $hic_matrices -n $sample_labels -chr $chrom -s $start_bin -e $stop_bin -r $bin_size -o $region -hist $bedgraph_files -hl $bedgraph_files_labels -hc FF0000,0000FF FF0000,0000FF -fh $fileheader -ext $filetype -high $highlight_opt -c $compare -p $pair -ptr 1 -trh 25 -spi 1 -si 1 -pdb $domainbars -pcd 1 -pcdf $domain_files #-g $gene_path
    python $hicplotter_abs_path -v -sif 1 -f $hic_matrices -n $sample_labels -chr $chrom -s $start_bin -e $stop_bin -r $bin_size -o $region -hist $bedgraph_files -hl $bedgraph_files_labels -hc FF0000,0000FF FF0000,0000FF -fh $fileheader -ext $filetype -high $highlight_opt -c $compare -p $pair -ptr 1 -trh 25 -spi 1 -si 1 -pdb $domainbars -pcd 1 -pcdf $domain_files #-g $gene_path
    foreach fout (chr*.$filetype)
      mv -f $fout $p/$outdir/`echo $fout | tr ':' ' ' | tr '-' ' ' | awk -F" " '{print $1"_"$2"-"$3}'`.$filetype
    end
  end
  cd $p
end

# Cleanup
rm -rf $workdir

goto_done:

# -------------------------------------
# -----  MAIN CODE ABOVE --------------
# -------------------------------------


# save variables
source ./code/code.main/scripts-save-job-vars

# done
scripts-send2err "Done."


