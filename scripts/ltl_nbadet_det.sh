#!/bin/sh
cat $1 | ltl2tgba -B | ../build/bin/nbadet -d -a -b -c -n -s
