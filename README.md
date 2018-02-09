# nbautils

## Getting the dependencies

First initialize and download the dependencies:
```
git submodule update --init --recursive --remote

wget http://automata.tools/hoa/cpphoafparser/down/cpphoafparser-0.99.2.tgz
tar xf cpphoafparser-0.99.2.tgz -C vendor/
mv vendor/cpphoafparser{-0.99.2,}
```

## Building nbautils

When all dependencies are satisfied, building is easy:

```
cmake CMakeLists.txt
make
```

## Contributing

Please run `make clangformat` before pushing code or issuing a pull request.

