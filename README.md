# nbautils

## Getting the dependencies

First initialize and download the dependencies:
```
git submodule update --init --recursive --remote
cd vendor/range-v3
git checkout 0.4.0
cd ../..
```

## Building nbautils

When all dependencies are satisfied, building is easy:

```
cmake CMakeLists.txt
make
```

## Using nbadet

The nbadet tool can be used like any other unix-style tool by piping the input automaton
into it. The current implementation only supports state-based Büchi automata as input, but
provides transition-based parity automata as output. Automata can be easily converted
using e.g. the `autfilt` tool from the [spot](https://spot.lrde.epita.fr/) library.
Without any flags, `nbadet` performs a completely unoptimized determinization. The
following example illustrates usage with a reasonable set of optimizations:

```
ltl2tgba -B 'FG b | (FG c & GF a)' | nbadet -k -j -t -i -r -o -m -d -u2
```

## Contributing

Please run `make clangformat` before pushing code or issuing a pull request.

