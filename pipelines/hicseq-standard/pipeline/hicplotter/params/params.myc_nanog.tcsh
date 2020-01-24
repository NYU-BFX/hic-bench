#!/bin/tcsh

source ./inputs/params/params.tcsh

module unload gcc               # this is necessary in order to take care of module conflicts in our system
module unload python
module load python/cpu/2.7.15-ES

# HiCplotter path
set hicplotter_path = ./code/HiCPlotter2.py
set hicplotter_params = "-hmc 1 -ptr 1"

# create bedgraphs for boundary scores
set branch_short = `echo $branch | sed 's/.*results\///' | sed 's/^matrix-distnorm.[^/]\+\///'`
set group_var = `echo $branch_short | cut -d'/' -f1 | cut -d'.' -f2`
set bscores_branch = ../boundary-scores/results/boundary-scores.$group_var.activity_500kb/$branch_short
set cell_type = `./code/code.main/read-sample-sheet.tcsh $sheet "$objects" cell-type`
set methods = (inter DI ratio)
set kappas = `cd $bscores_branch/$objects[1]; ls -1 all_scores.k=*.tsv | cut -d'.' -f2`
set bedgraphs = ()
set bedgraph_labels = ()
foreach kappa ($kappas)
  set f = $bscores_branch/$objects[1]/all_scores.$kappa.tsv
  set bfiles = ()
  set blabels = ()
  
  # add bscores for each method
  foreach m ($methods)
    set k = `head -1 $f | tr '\t' '\n' | grep -n "^$m"'$' | cut -d':' -f1`             # get column number for method
    cat $f | sed '1d' | cut -f1,$k | sed 's/:/\t/' | sed 's/-/\t/' | grep -v '	NA$' >! $outdir/bscores.$kappa.$m.bedGraph               # convert to bedgraph, remove NA values
    set bfiles = ($bfiles ../bscores.$kappa.$m.bedGraph)
    set blabels = ($blabels $m)
  end

  # add CTCF ChIP-seq
  if (-e inputs/data.external/$cell_type/CTCF.bedGraph) then
    set bfiles = ($bfiles `readlink -f inputs/data.external/$cell_type/CTCF.bedGraph`)
    set blabels = ($blabels CTCF)
  endif

  # add to the lists
  set bedgraphs = ($bedgraphs `echo $bfiles | tr ' ' ','`)
  set bedgraph_labels = ($bedgraph_labels `echo $blabels | tr ' ' ','`)
end

# regions to plot
set genes = 'MYC|NANOG'
set regions = `cat $genome_dir/gene-name.bed | grep -wiE "$genes" | gtools-regions center | gtools-regions shiftp -5p -4000000 -3p +4000000 | cut -f-3 | sed 's/\t/:/' | sed 's/\t/-/'`
set tiles = "$outdir/regions.bed"
cat $genome_dir/gene-name.bed | grep -wiE "$genes" | sed 's/^/0.7\t66,80,209\t/' | tools-cols -t 2 3 4 0 1 5 >! $tiles
set tiles_labels = "regions"
set highlight = 0
set highlight_bed = ""
set fileheader = 0         # Either 1 or 0 (header / no header)
set insulation_score = 0   # Either 1 or 0 (include insulation index or not)



