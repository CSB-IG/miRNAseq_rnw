#!/bin/bash

# prints pruned to STDOUT
# usage: pruner.sh THRESHOLD PRUNEE
# prune by mi value, ej. 0.101
# for table: node node interaction ($3 = mi value)
THRESHOLD=$1

awk -F'\t' "\$3 > $THRESHOLD {print \$1\"\t\"\$2\"\t\"\$3}" $2