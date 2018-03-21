#!/bin/sh
pat="aut_XXXX.hoa"
curr=$(mktemp /tmp/$pat)
lines=0
while read line; do
  if [ $lines == "0" ]; then
    echo $curr
  fi

  echo $line >> $curr
  lines=$(($lines+1))

  if [ "$line" == "--END--" ]; then
    curr=$(mktemp /tmp/$pat)
    lines=0
  fi

done
