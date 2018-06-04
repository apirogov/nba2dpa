#!/bin/sh

for SZ in '5..10' '10..20' '20..30' '30..50' '50..75' '75..100'; do
  for ACP in 0.2 0.4 0.6 0.8; do
    for DEN in 0.2 0.4 0.6 0.8; do
      randaut -u -B -a $ACP -e $DEN -n -1 -Q $SZ 2..6 |
        autfilt --is-deterministic -v -n 100
    done
  done
done
