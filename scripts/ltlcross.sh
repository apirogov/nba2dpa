#!/bin/sh
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
wrap() {
  echo 'cat %F |' $SCRIPTPATH/$1 '> %O'
}
echo $(wrap ./ltl_dstar.sh)
ltlcross -T 180 --strength \
  -t "$(wrap ltl_dstar.sh)" \
  -t "$(wrap ltl_spot.sh)" \
  -t "$(wrap ltl_rabinizer.sh)" \
  -t "$(wrap ltl_nbadet.sh)" \
 --csv=$1
  #-t "cat %F | $ltl_nbadet | autfilt --small -S -C > %O" \
