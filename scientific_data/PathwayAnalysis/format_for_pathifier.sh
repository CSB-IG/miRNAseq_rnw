#!/bin/bash
# turns a gmt into a squared matrix for pathifier adding tabs into blank spaces

file="$1"
cols=$(awk '{if(max<NF) {max=NF}} END {print max}' "$file")
awk '{ for( i=1; i<391; i++) {printf( $i "\t")}; printf($'"$cols"' "\n")} ' "$file" 
