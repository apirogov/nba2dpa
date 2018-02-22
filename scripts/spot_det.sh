#!/bin/sh
cat $1 | autfilt --complement | autfilt --complement --small --high -S
