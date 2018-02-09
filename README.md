# nbautils

## Getting the dependencies

First initialize and download the dependencies:
```
git submodule update --init --recursive --remote
```

## Building nbautils

When all dependencies are satisfied, building is easy:

```
cmake CMakeLists.txt
make
```

## Contributing

Please run `make clangformat` before pushing code or issuing a pull request.

