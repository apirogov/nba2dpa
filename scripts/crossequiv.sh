#!/bin/bash
SEED=$(date +%s)
randltl --seed $SEED --tree-size 10..25 4 -n 100 | while read ltl; do
  echo $ltl
  SRC=$(mktemp)
  echo $ltl | ltl2tgba -B > $SRC

  FILES="$SRC"
  for cmd in "$@"; do
    RES=$(mktemp)
    cat $SRC | $cmd > $RES
    FILES="$FILES $RES"
  done
  # echo $FILES
  for a in $(echo $FILES); do
    # echo a: $a
    for b in $(echo $FILES); do
      # echo b: $b
      echo -n "$(autfilt -c --equivalent-to=$(realpath $b) $a) "
    done
  done
  echo ""
done
