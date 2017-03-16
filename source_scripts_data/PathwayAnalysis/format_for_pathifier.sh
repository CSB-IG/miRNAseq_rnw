#!/bin/bash
# turns a gmt file into a squared matrix for pathifier adding tabs into blank spaces
# in order to run the pathifier R script correctly

file="$1"
cols=$(awk '{if(max<NF) {max=NF}} END {print max}' "$file")
awk '{ for( i=1; i<391; i++) {printf( $i "\t")}; printf($'"$cols"' "\n")} ' "$file" 
