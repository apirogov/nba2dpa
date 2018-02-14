cmake_minimum_required(VERSION 3.1)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)

project(nbautils)

set(nbautils-build_tests OFF CACHE BOOL "Whether to build tests")

# Enable C++14
set (CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE Release)
#endif()


##################################
# Find thread
find_package(Threads REQUIRED)

# Find boost
#find_package(Boost 1.40 REQUIRED)
#include_directories(${Boost_INCLUDE_DIR})

# Catch for testing
set(Catch_INCLUDE_DIR vendor/Catch/single_include)

# Spdlog for logging
set(spdlog_SOURCE_DIR vendor/spdlog)
include_directories(${spdlog_SOURCE_DIR}/include)

# args for arg parsing
set(args_SOURCE_DIR vendor/args)
include_directories(${args_SOURCE_DIR})

# cpphoafparser
#the following is used to prevent re-downloading on clean build, adapted from:
#https://github.com/anton-sergeev/cmake_externalproject/blob/b3ef26fe616ac98a3b3eeb86a6b33add0765cacf/CMakeLists.txt
set(EXTERNALS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor)
set(CPPHOAF_VERSION "0.99.2")
set(CPPHOAF_TARBALL_URL http://automata.tools/hoa/cpphoafparser/down/cpphoafparser-${CPPHOAF_VERSION}.tgz)
set(CPPHOAF_TARBALL ${CMAKE_CURRENT_SOURCE_DIR}/tarballs/cpphoafparser-${CPPHOAF_VERSION}.tgz)
set(CPPHOAF_DIR ${EXTERNALS_DIR}/cpphoafparser)
if(NOT EXISTS "${CPPHOAF_TARBALL}")
  file(DOWNLOAD ${CPPHOAF_TARBALL_URL} ${CPPHOAF_TARBALL})
endif()
add_custom_target(cpphoaf-download DEPENDS ${CPPHOAF_TARBALL})

include(ExternalProject)
ExternalProject_Add(cpphoafparser
DEPENDS cpphoaf-download
URL ${CPPHOAF_TARBALL}
URL_MD5 "9a36d26d42603782bdfbcd01c081bc75"
CMAKE_ARGS ""
PREFIX ${CPPHOAF_DIR}
UPDATE_COMMAND cmake -E copy_directory ${EXTERNALS_DIR}/cpphoafparser-${CPPHOAF_VERSION}-overlay/. ${CPPHOAF_DIR}
CONFIGURE_COMMAND ""
BUILD_COMMAND ""
INSTALL_COMMAND ""
#TEST_COMMAND ""
)
ExternalProject_Get_Property(cpphoafparser source_dir)
include_directories(${source_dir}/include)
##################################

include_directories(src)

set(nbautils_SOURCE src/memusage.h src/memusage.c
                   src/bench.hh src/bench.cc
                   src/debug.hh src/debug.cc
                   src/relorder.hh src/relorder.cc
                   src/types.hh src/triebimap.hh
                   src/io.hh src/io.cc
                   src/ps.hh #src/ps.cc
                   src/scc.hh #src/scc.cc
                   src/level.hh src/level.cc
                   src/det.hh src/det.cc
                   )

add_library(nbautils-lib STATIC ${nbautils_SOURCE})
target_link_libraries(nbautils-lib ${CMAKE_THREAD_LIBS_INIT})
add_dependencies(nbautils-lib cpphoafparser)

add_executable(scratchpad EXCLUDE_FROM_ALL src/scratchpad.cc)
target_link_libraries(scratchpad nbautils-lib ${CMAKE_THREAD_LIBS_INIT})

add_executable(nbadet src/nbadet.cc)
target_link_libraries(nbadet nbautils-lib ${CMAKE_THREAD_LIBS_INIT})

if (nbautils-build_tests)
    enable_testing()

    add_executable(nbautils-test test/nbautils_test.cc test/test_main.cc)
    target_include_directories(nbautils-test PUBLIC ${Catch_INCLUDE_DIR})
    target_link_libraries(nbautils-test nbautils-lib)

    # to add binaries to "make test" suite
    macro (create_test target)
      add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}> WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test)
    endmacro (create_test)

    create_test(nbautils-test)
endif()

add_custom_target(
        clangformat
        COMMAND /usr/bin/clang-format
        -style=file
        -i
        src/*
)

