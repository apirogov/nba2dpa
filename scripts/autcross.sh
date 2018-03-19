#!/bin/sh
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
wrap() {
  echo 'cat %H |' $@ '> %O'
}
randaut -B -n 100 1..5 | autcross --language-preserved -T 180 \
  -t "$(wrap ltl2dstar -B -H - -)" \
  -t "$(wrap autfilt -D -S -C -G)" \
  -t "$(wrap $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -c -t -p -m)" \
 $@
