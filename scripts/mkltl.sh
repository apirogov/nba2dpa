#!/bin/bash
SEED=$(date +%s)
WANT=$3
NUM=0
while [ $NUM -lt $WANT ]; do
  BATCH=$(($WANT-$NUM))
  if [ $BATCH -gt 10 ]; then
    BATCH=10
  fi

  randltl --seed $SEED --tree-size $2 -n 100 $1 | awk '{ print length(), $0 | "sort -n" }' |
    cut -d' ' -f 2- | tail -n $BATCH

  SEED=$(($SEED+1))
  NUM=$(($NUM+$BATCH))
done
