# nbautils

## Getting the dependencies

First initialize the directories of dependencies:
```
git submodule update --init --recursive --remote
```

To build spot, ensure you have flex and bison installed and a modern g++/clang toolchain.
```
cd vendor/spot
autoreconf -vfi
./configure --disable-python --disable-shared --enable-static --prefix=$(realpath ..)/spot-build
make
make install
```

## Building nbautils

When all dependencies are satisfied, building is easy:

```
cmake CMakeLists.txt
make
```

## Contributing

Please run `make clangformat` before pushing code or issuing a pull request.

