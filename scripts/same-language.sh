#!/bin/bash
#pipe automaton into tool stdin, check language equivalence with result on stdout
#usage: ./same-language.sh someaut.hoa usual command with parameters
# FILE=$(realpath $1)
FILE=$(mktemp /dev/shm/XXXX.hoa)
cat > $FILE
cat $FILE | $@ | autfilt -c --equivalent-to=<(cat $FILE | autfilt --trust-hoa=no -D -P)
