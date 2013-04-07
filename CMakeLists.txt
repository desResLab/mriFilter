project(mpFilterApp)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
include_directories(".")
add_definitions("-Wall -m64 -O3 -std=c++0x")
add_executable(${PROJECT_NAME} ${SRC_LIST})





