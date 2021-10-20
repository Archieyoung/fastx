cmake_minimum_required(VERSION 3.10)
# project and version
project(GLAT LANGUAGES CXX C VERSION 0.0.1.0)

include_directories(vendor/htslib)

# PROJECT_BINARY_DIR default is the build dir
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_BUILD_TYPE "Release")

# CXX configure
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    message(STATUS "[STATUS]: Build Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -pedantic -fprofile-arcs -ftest-coverage")
elseif (CMAKE_BUILD_TYPE STREQUAL Release)
    message(STATUS "[STATUS]: Build Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

find_package(Boost 1.74 REQUIRED COMPONENTS filesystem)
include_directories(${Boost_INCLUDE_DIRS})

# htslib library and include files in vendor
set(HTS_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
set(HTS_LIB ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib/libhts.a)

if (NOT EXISTS ${HTS_INCLUDE})
    message(FATAL_ERROR "[FATAL_ERROR]: Can not find htslib directory in \
${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib, Please check vendor directory \
for existence of htslib.")
endif (NOT EXISTS ${HTS_INCLUDE})

if (EXISTS ${HTS_LIB})
    message(STATUS "[STATUS]: Find htslib static library libhts.a")
else (EXISTS ${HTS_LIB})
    message(STATUS "[STATUS]: htslib static library libhts.a not exist in \
${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib, build one ...")
    execute_process(COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
endif (EXISTS ${HTS_LIB})


add_executable(fastx
    src/common.cpp
    src/fastx_head.cpp
    src/fastx_sample.cpp
    src/fastx.cpp)

target_include_directories(fastx PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib
    $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS})

target_link_libraries(fastx
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib/libhts.a
    z m bz2 lzma pthread curl ${Boost_LIBRARIES})