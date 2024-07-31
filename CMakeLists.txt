cmake_minimum_required(VERSION 3.29)
project(openkraken)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/)
include_directories(${CMAKE_SOURCE_DIR}/../)

file(GLOB_RECURSE openkraken_sources
    "./openkraken/" "*.cpp" "*.hpp"
)

add_definitions(-DINTERFACEDLL_EXPORTS)
add_library(openkraken SHARED ${openkraken_sources})
