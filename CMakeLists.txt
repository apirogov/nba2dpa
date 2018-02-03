cmake_minimum_required(VERSION 3.1)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)

project(nbautils)

set(nbautils-build_tests OFF CACHE BOOL "Whether to build tests")

# Enable C++14
set (CMAKE_CXX_STANDARD 14)

# Find thread
find_package(Threads REQUIRED)

# Catch for testing
set(Catch_INCLUDE_DIR vendor/Catch/single_include)

# Spdlog for logging
set(spdlog_SOURCE_DIR vendor/spdlog)
include_directories(${spdlog_SOURCE_DIR}/include)

# args for arg parsing
set(args_SOURCE_DIR vendor/args)
include_directories(${args_SOURCE_DIR})

##################################
# spot
set(spot_DIST_DIR vendor/spot-build)
include_directories(${spot_DIST_DIR}/include)
add_library(spot STATIC IMPORTED) #add imported target
set_property(TARGET spot
             PROPERTY IMPORTED_LOCATION ${spot_DIST_DIR}/lib/libspot.a)
# spot brings a version of bddx it requires
add_library(bddx STATIC IMPORTED)
set_property(TARGET bddx
             PROPERTY IMPORTED_LOCATION ${spot_DIST_DIR}/lib/libbddx.a)

##################################
add_custom_target(
        clangformat
        COMMAND /usr/bin/clang-format
        -style=file
        -i
        src/*
)

include_directories(src)

set(nbautils_SOURCE src/relorder.hh src/relorder.cc
                   src/nba.hh src/nba.cc
                   src/dpa.hh src/dpa.cc
                   src/level.hh)

add_library(nbautils-lib STATIC ${nbautils_SOURCE})
#add_executable(nbautils ${nbautils_SOURCE})
target_link_libraries(nbautils-lib ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(nbautils-lib spot)

add_executable(scratchpad src/scratchpad.cc ${nbautils_SOURCE})
target_link_libraries(scratchpad ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(scratchpad spot bddx)

if (nbautils-build_tests)
    enable_testing()

    add_executable(nbautils-test test/nbautils_test.cc test/test_main.cc)
    target_include_directories(nbautils-test PUBLIC ${Catch_INCLUDE_DIR})
    target_link_libraries(nbautils-test nbautils-lib)
    #target_link_libraries(nbautils-test spot)

    add_test(nbautils_test nbautils-test)
endif()
