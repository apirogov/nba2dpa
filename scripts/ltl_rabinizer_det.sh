#!/bin/sh
cat $1 | ltl2dpa --mode=rabinizer | autfilt -S
