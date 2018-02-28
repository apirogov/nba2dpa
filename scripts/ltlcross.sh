#!/bin/sh
spot_ltl2dpa="ltl2tgba | autfilt --complement | autfilt --complement --small --high -S"
ltl_nbadet="ltl2tgba -B | ../build/bin/nbadet -d -s -n -a -b -c"
ltlcross \
  -t "cat %F | $spot_ltl2dpa > %O" \
  -t "cat %F | $ltl_nbadet   > %O" \
  # -t "cat %F | ltl2tgba -B | ../build/bin/nbadet    > %O" \
#  --csv=$2

  #-t "cat %F | ./ltl_nbadet_det.sh | autfilt --small -S   > %O" \
  #-t "cat %F | ./ltl_dstar_det.sh     > %O" \
  #-t "cat %F | ./ltl_rabinizer_det.sh > %O"
