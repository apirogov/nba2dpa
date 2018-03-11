#!/bin/sh
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
ltl2tgba -B | $SCRIPTPATH/../build/bin/nbadet -d -s -n -a -b -c -t -p -m $@  #-d -s -n -a -b -c -t -p -m
