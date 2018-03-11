#!/bin/sh
cat $1 | ltl2tgba -B | autfilt --complement | autfilt --complement --small -S -C
