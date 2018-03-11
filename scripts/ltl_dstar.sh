#!/bin/sh
cat $1 | ltl2tgba -B | ./dstar_det.sh # | autfilt -S -C
