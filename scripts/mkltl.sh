#/bin/sh
randltl --seed $(date +%s) --tree-size $2 -n 100 $1 | awk '{ print length(), $0 | "sort -n" }' | cut -d' ' -f 2- | tail -n $3
