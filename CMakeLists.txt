cmake_minimum_required(VERSION 3.10)
# project and version
project(GLAT LANGUAGES CXX C VERSION 0.0.1.0)

include_directories(vendor/htslib)

# PROJECT_BINARY_DIR default is the build dir
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (NOT EXISTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    execute_process(COMMAND mkdir -p ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif (NOT EXISTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if (NOT EXISTS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    execute_process(COMMAND mkdir -p ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif (NOT EXISTS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# CXX configure
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    message(STATUS "[STATUS]: Build Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++17 -g -Wall -Wextra -pedantic -fprofile-arcs -ftest-coverage")
elseif (CMAKE_BUILD_TYPE STREQUAL Release)
    message(STATUS "[STATUS]: Build Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++17 -O3 -Wall -Wextra -pedantic")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

set(MAKE_JOBS 4 CACHE STRING "number of jobs for run make in this cmake list")

set(USE_LIBDEFLATE ON CACHE STRING "if use libdeflate in htslib?")

# libdeflate library and include files in vendor
set(LIBDEFLATE ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdeflate)
set(LIBDEFLATE_LIB ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdeflate/libdeflate.so)

if (${USE_LIBDEFLATE})
    if (NOT EXISTS ${LIBDEFLATE})
        message(FATAL_ERROR "[FATAL_ERROR]: Can not find libdeflate directory in \
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor, Please check vendor directory \
    for existence of libdeflate.")
    endif (NOT EXISTS ${LIBDEFLATE})
endif (${USE_LIBDEFLATE})

if (${USE_LIBDEFLATE})
    if (EXISTS ${LIBDEFLATE_LIB})
        message(STATUS "[STATUS]: Find libdeflate library libdeflate.so")
    else (EXISTS ${LIBDEFLATE_LIB})
        message(STATUS "[STATUS]: libdeflate library libdeflate.so not exist in \
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdeflate, build one ...")
        execute_process(COMMAND make -j ${MAKE_JOBS} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdeflate)
    endif (EXISTS ${LIBDEFLATE_LIB})
    execute_process(COMMAND cp ${LIBDEFLATE}/libdeflate.a ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    execute_process(COMMAND cp ${LIBDEFLATE}/libdeflate.so ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    execute_process(COMMAND cp ${LIBDEFLATE}/libdeflate.so.0 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif (${USE_LIBDEFLATE})


# htslib library and include files in vendor
set(HTSLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
set(HTSLIB_LIB ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib/libhts.so)

if (NOT EXISTS ${HTSLIB_DIR})
        message(FATAL_ERROR "[FATAL_ERROR]: Can not find htslib directory in \
${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib, Please check vendor directory \
for existence of htslib.")
endif (NOT EXISTS ${HTSLIB_DIR})


if (EXISTS ${HTSLIB_LIB})
    message(STATUS "[STATUS]: Find htslib shared library libhts.so")
else (EXISTS ${HTSLIB_LIB})
    message(STATUS "[STATUS]: htslib shared library libhts.so not exist in ${HTSLIB_DIR}, build one ...")
    execute_process(COMMAND autoheader WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    execute_process(COMMAND autoconf WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    if (${USE_LIBDEFLATE})
        execute_process(COMMAND bash -c "./configure CPPFLAGS=-I${LIBDEFLATE} LDFLAGS='-L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -Wl,-R${CMAKE_LIBRARY_OUTPUT_DIRECTORY}' --with-libdeflate" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    else (${USE_LIBDEFLATE})
        execute_process(COMMAND bash -c "./configure" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    endif (${USE_LIBDEFLATE})
    execute_process(COMMAND make -j ${MAKE_JOBS} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    execute_process(COMMAND cp ${HTSLIB_LIB} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    execute_process(COMMAND cp ${HTSLIB_LIB}.3 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif (EXISTS ${HTSLIB_LIB})

set(ZLIB_COMPAT ON)
add_subdirectory(vendor/zlib-ng)

add_executable(fastx
    src/common.cpp
    src/fastx_head.cpp
    src/fastx_sample.cpp
    src/fastx_split.cpp
    src/fastx.cpp)

target_include_directories(fastx PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib
    )

target_link_libraries(fastx
    ${HTSLIB_LIB}
    zlibstatic m bz2 lzma pthread curl)