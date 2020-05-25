#!/bin/tcsh

source ./inputs/params/params.tcsh

set APA_resolution = 10000		     # set APA resolution
set URm = 2.5			             # Upper Right Multiplier: to set the maximum value in the heatmap color scale [ UR = URm * mean(upper right quadrant counts) ]    
