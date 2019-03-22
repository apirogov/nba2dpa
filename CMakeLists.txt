cmake_minimum_required(VERSION 3.1)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/bin)

project(nbautils)

set(nbautils-build_tests OFF CACHE BOOL "Whether to build tests")

# Enable C++14
set (CMAKE_CXX_STANDARD 14)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fconcepts HAS_FCONCEPTS)
if(HAS_FCONCEPTS)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -fconcepts")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")
endif()

#cmake -DCMAKE_BUILD_TYPE=Release|Debug
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

# rapidcheck for testing
#add_subdirectory("vendor/rapidcheck")

# range-v3
set(rangev3_INCLUDE_DIR vendor/range-v3/include)
include_directories(${rangev3_INCLUDE_DIR})

# Spdlog for logging
set(spdlog_SOURCE_DIR vendor/spdlog)
include_directories(${spdlog_SOURCE_DIR}/include)

# args for arg parsing
set(args_SOURCE_DIR vendor/args)
include_directories(${args_SOURCE_DIR})

# process lib
set(procxx_SOURCE_DIR vendor/procxx)
include_directories(${procxx_SOURCE_DIR}/include)

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

set(nbautils_SOURCE
                   src/metrics/memusage.h src/metrics/memusage.c src/metrics/bench.hh
                   src/common/bimap.hh src/common/util.hh src/common/scc.hh
                   src/common/types.hh src/common/types.cc
                   src/io.hh src/io.cc
                   src/aut.hh src/ps.hh
                   src/det.hh src/det.cc
                   src/detstate.hh src/detstate.cc
                   src/pa.hh src/pa.cc

				   # Complementation of Buechi-automata #################
					src/compl/compl_tag.hh
					src/compl/compl_print.hh
					src/compl/compl_constr1.hh
					src/compl/compl_constr2.hh
					src/compl/compl_constr3.hh
				   ######################################################
                   )

add_library(nbautils-lib STATIC ${nbautils_SOURCE})
target_link_libraries(nbautils-lib ${CMAKE_THREAD_LIBS_INIT})
add_dependencies(nbautils-lib cpphoafparser)

add_executable(nbadet src/tools/nbadet.cc)
target_link_libraries(nbadet nbautils-lib ${CMAKE_THREAD_LIBS_INIT})

#add_executable(nba2dot EXCLUDE_FROM_ALL src/nba2dot.cc)
#target_link_libraries(nba2dot nbautils-lib ${CMAKE_THREAD_LIBS_INIT})

#add_executable(nbasccinfo EXCLUDE_FROM_ALL src/nbasccinfo.cc)
#target_link_libraries(nbasccinfo nbautils-lib ${CMAKE_THREAD_LIBS_INIT})

add_executable(permsets EXCLUDE_FROM_ALL src/tools/permsets.cc)
target_link_libraries(permsets nbautils-lib ${CMAKE_THREAD_LIBS_INIT})

add_executable(scratchpad EXCLUDE_FROM_ALL src/tools/scratchpad.cc)
target_link_libraries(scratchpad nbautils-lib ${CMAKE_THREAD_LIBS_INIT})


# Complementation of Buechi-automata ################################################
add_executable(compl EXCLUDE_FROM_ALL src/tools/compl.cc)							#
target_link_libraries(compl nbautils-lib ${CMAKE_THREAD_LIBS_INIT})					#
#####################################################################################

if (nbautils-build_tests)
    enable_testing()
    # to add binaries to "make test" suite
    macro (create_test target)
      add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
    endmacro (create_test)

    add_library(test-main test/test_main.cc)
    target_include_directories(test-main PUBLIC ${Catch_INCLUDE_DIR})
    target_link_libraries(test-main nbautils-lib) #rapidcheck)

    set(test_nbautils_SOURCE
      #                       test/test_nbautils_common.cc
      #                       test/test_nbautils_swa.cc
      #                       test/test_nbautils_scc.cc
      #                       test/test_nbautils_ps.cc
      #                       test/test_nbautils_pa.cc
                            )

    # add_executable(test-nbautils ${test_nbautils_SOURCE})
    # target_include_directories(test-nbautils PUBLIC ${Catch_INCLUDE_DIR}) # ${rapidcheck_INCLUDE_DIR})
    # target_link_libraries(test-nbautils test-main)

    # create_test(test-nbautils)
endif()

add_custom_target(
        clangformat
        COMMAND /usr/bin/clang-format
        -style=file
        -i
        src/*
)

