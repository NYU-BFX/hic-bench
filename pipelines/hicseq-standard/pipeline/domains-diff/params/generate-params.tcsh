#!/bin/tcsh

foreach tad_caller (hicratio.d_0500 topdom.win_5)
  foreach norm (cpm dist_norm)
    foreach ref1 (FALSE TRUE)
      if ($ref1 == FALSE) then
        set ref1_name = common
      else
        set ref1_name = ref1
      endif

      set p = params.$tad_caller.$norm.$ref1_name.tcsh
      
      echo "Generating parameter script $p..."

      echo '#\!/bin/tcsh' >! $p
      echo '' >> $p
      echo "source ./inputs/params/params.tcsh" >> $p
      echo '' >> $p
      echo "set tad_caller = $tad_caller" >> $p
      echo "set is_normalize = $norm" >> $p
      echo "set use_sample1_ref = $ref1" >> $p
      echo '' >> $p
      echo 'source params/params-template.tcsh' >> $p

      chmod +x $p
    end
  end
end



