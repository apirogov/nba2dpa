#/bin/sh
randltl --seed $(date +%s) --tree-size 50 5 | ltl2tgba -B | tee $1
