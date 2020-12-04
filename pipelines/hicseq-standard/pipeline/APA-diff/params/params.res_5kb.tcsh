#!/bin/tcsh

source ./inputs/params/params.tcsh

set APA_resolution = 5000		     # set APA resolution
set URm = 2.5			             # Upper Right Multiplier: to set the maximum value in the heatmap color scale [ UR = URm * mean(upper right quadrant counts) ]    
set fmin_distance = 12			     # minimimum distance factor (mdf): mdf * (APA_resolution/sqrt(2)) = minimum distance from diagonal (e.g: mdf=12; res=10000 -> min.distance=84853) 
