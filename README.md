# nbautils

This is the prototype implementation of the new determinization construction and
optimizations described in 
[ICALP'19](http://dx.doi.org/10.4230/LIPIcs.ICALP.2019.120) and 
[ATVA'19](https://doi.org/10.1007/978-3-030-31784-3_18).

## Getting the dependencies

First initialize and download the dependencies:
```
git submodule update --init --recursive --remote
cd vendor/range-v3
git checkout 0.4.0
cd ../..
cd vendor/spdlog
git checkout v0.16.3
ch ../..
```

## Building nbautils

When all dependencies are satisfied, building is easy (tested with CMake 3.19 and g++ 10.2):

```
cmake CMakeLists.txt
make
```

## Using nbadet

The nbadet tool can be used like any other unix-style tool by piping the input automaton
into it. The current implementation only supports state-based Büchi automata as input, but
provides transition-based parity automata as output. Automata can be easily converted
using e.g. the `autfilt` tool from the [spot](https://spot.lrde.epita.fr/) library.
Without any flags, `nbadet` performs a completely unoptimized determinization. 
A reasonable default set of optimizations is: `-k -j -t -i -r -o -m -d -u1`.

The following example illustrates usage:

```
genltl --ms-phi-h 3 | ltl2tgba -B | nbadet -j -k -t -i -r -a -b -m
```

## Contributing

Please run `make clangformat` before pushing code or issuing a pull request.

