#!/bin/sh
spot_ltl2dpa="ltl2tgba -B | autfilt --complement | autfilt --complement --small -S"
ltl_nbadet="ltl2tgba -B | ../build/bin/nbadet -d -s -n -a -b -t -p -m "

ltlcross \
  -t "cat %F | $spot_ltl2dpa > %O" \
  -t "cat %F | $ltl_nbadet   > %O" \
 --csv=$1
#  -t "cat %F | $ltl_nbadet | autfilt --small -S   > %O" \
  #-t "cat %F | ./ltl_dstar_det.sh     > %O" \
  #-t "cat %F | ./ltl_rabinizer_det.sh > %O"
