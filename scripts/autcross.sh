#!/bin/sh
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
wrap() {
  echo 'cat %H |' $@ '> %O'
}
# randaut --seed=1234 -u -B -n 500 1..3 -Q 4..10 |
cat | autcross -T 180 \
  -t "$(wrap autfilt --trust-hoa=no -D -S -C -p)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d )" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s -n -a)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -c)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -c -t)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -c -t -g)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -c -t -g -u1)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -c -t -g -u1 -l)" \
 $@
#  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -c -t -g -u1 -l -p -m)" \
  # -t "$(wrap ltl2dstar -B -H - -)" \
