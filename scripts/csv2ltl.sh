#!/bin/sh
#extract formulas back from CSV output
cat $1 | cut -d \, -f 1 | uniq | awk '{if(NR%2==0){print}}' | sed 's/"//g'
