cmake_minimum_required(VERSION 3.2)
project(spoa)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -march=native")

option(spoa_build_tests "Build spoa google tests" OFF)

add_executable(spoa
    src/alignment_engine.cpp
    src/chain.cpp
    src/graph.cpp
    src/main.cpp
    src/simd_alignment_engine.cpp
    src/sisd_alignment_engine.cpp)

target_include_directories(spoa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

add_subdirectory(vendor/bioparser)
target_link_libraries(spoa bioparser)

if (spoa_build_tests)
endif(spoa_build_tests)

install(TARGETS spoa DESTINATION bin)
