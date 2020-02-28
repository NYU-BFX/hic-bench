#!/bin/tcsh

source ./inputs/params/params.tcsh

set format = juicer

module unload r                  # there is a conflict with juicer on our system...
module load juicer/1.5

