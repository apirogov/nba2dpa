#!/bin/bash
#pipe automaton into tool stdin, check language equivalence with result on stdout
#usage: ./same-language.sh someaut.hoa usual command with parameters
FILE=$(realpath $1)
cat $FILE | ${@:2} | autfilt -c --equivalent-to=$FILE
