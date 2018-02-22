#!/bin/sh
ltlcross -F $1 \
  -t "cat %F | ./ltl_spot_det.sh      > %O" \
  -t "cat %F | ./ltl_nbadet_det.sh    > %O" \
  -t "cat %F | ./ltl_nbadet_det.sh | autfilt --small    > %O" \
  --csv=$2
  #-t "cat %F | ./ltl_dstar_det.sh     > %O" \
  #-t "cat %F | ./ltl_rabinizer_det.sh > %O"
