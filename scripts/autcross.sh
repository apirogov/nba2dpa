#!/bin/bash
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
SEED=$(date +%s)
DATE=`date +%Y%m%d-%H%M%S`
wrap() {
  echo 'cat %H |' $@ '> %O'
}

# randaut --seed=1234 -u -B -n 500 1..3 -Q 4..10 |
# cat <(echo -e '0\n1') <(randltl --seed $SEED --tree-size 15..20 4 -n 1000 -r3) |
#   ltlfilt -u --remove-wm --size 40..70 | head -n 300 |
#   awk '{ print length($0) " " $0; }' $file | sort -n | cut -d ' ' -f 2- |
#   ltl2tgba -B |
cat |
  autcross \
  -t "$(wrap autfilt -D -C -P)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u0 -k -j -t -n -a)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u0 -k -j -t -n -a -b)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u0 -k -j -t -n -a -b -c)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u0 -k -j -t -n -a -b -c -e)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u0 -k -j -t -n -a -b -c -d -e)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u0 -k -j -t -n -a -b -c -d -e -l)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u1 -k -j -t -n -a)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u1 -k -j -t -n -a -b)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u1 -k -j -t -n -a -b -c)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u1 -k -j -t -n -a -b -c -e)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u1 -k -j -t -n -a -b -c -d -e)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u1 -k -j -t -n -a -b -c -d -e -l)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u2 -k -j -t -n -a)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u2 -k -j -t -n -a -b)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u2 -k -j -t -n -a -b -c)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u2 -k -j -t -n -a -b -c -e)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u2 -k -j -t -n -a -b -c -d -e)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -u2 -k -j -t -n -a -b -c -d -e -l)" \
  $@ --save-bogus=failed_$DATE.hoa --csv=stats_$DATE.csv
 # TODO: -t, -p -m
  # -t "$(wrap ltl2dstar -B -H - -)" \
