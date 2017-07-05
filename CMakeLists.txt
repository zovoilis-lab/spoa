cmake_minimum_required(VERSION 3.2)
project(spoa)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -march=native")

option(spoa_build_tests "Build spoa google tests" OFF)

add_library(spoa STATIC
    src/alignment_engine.cpp
    src/graph.cpp
    src/simd_alignment_engine.cpp
    src/sisd_alignment_engine.cpp)

add_executable(spoa-bin
    src/chain.cpp
    src/main.cpp)

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

target_link_libraries(spoa-bin spoa bioparser)
set_target_properties(spoa-bin PROPERTIES OUTPUT_NAME spoa)

target_include_directories(spoa PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

if (spoa_build_tests)
endif(spoa_build_tests)

install(TARGETS spoa-bin DESTINATION bin)
install(TARGETS spoa DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/spoa DESTINATION include)
