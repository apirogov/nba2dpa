#!/bin/sh
ltl2tgba | autfilt --complement | autfilt --complement --small --high -S
