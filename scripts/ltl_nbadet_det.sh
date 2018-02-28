#!/bin/sh
cat $1 | ltl2tgba -B | ../build/bin/nbadet -t -vv
