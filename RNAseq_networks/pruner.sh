
# prune by mi value, ej. 0.101
# for table: node node interaction ($3 = mi value)

awk -F'\t' '$3 > 0.101 {print $1"\t"$2"\t"$3}' listota.txt > listita.txt                    