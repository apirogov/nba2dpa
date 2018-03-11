#!/bin/bash
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
SEED=$(date +%s)
wrap() { echo 'cat %F |' $SCRIPTPATH/$1 '> %O'; }

#TODO: fix -b
#randltl --seed $SEED --tree-size 10..25 4 -n 100 |
cat <(echo -e '0\n1') <(randltl --seed $SEED --tree-size 15..20 4 -n 1000 -r3) | ltlfilt -u --remove-wm --size 1..70 | head -n 300 |
  stdbuf -i0 -e0 -o0 ltlcross -T 180 -t "ltl2tgba -D -G -S -C -F %F > %O" \
  -t "ltl2tgba -B -F %F | $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -c -t -p -m > %O" \
  --products=5 $@ 2>&1 |
  grep --line-buffered -v Running | grep --line-buffered -v Performing
