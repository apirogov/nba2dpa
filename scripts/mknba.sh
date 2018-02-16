#!/bin/sh
SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
wantstates=$1
filename=$2

twicewant=$(($wantstates * 2))
while true; do
  $SCRIPTPATH/mkltl.sh 5 20 100 | while read formula; do
    echo $formula | ltl2tgba -B > $filename
    numst=$(grep States $filename | awk '{print $2}')
    # echo $numst
    # TODO: ignore deterministic
    if [ "$numst" -ge "$wantstates" -a "$numst" -le "$twicewant" ]; then
      exit 1
    fi
  done
  if [ "$?" = "1" ]; then
    break
  fi
done
cat $filename
