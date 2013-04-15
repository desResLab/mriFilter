project(mpFilterApp)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
include_directories(".")

message("Including Boost Folders...")
INCLUDE_DIRECTORIES("G:/Daniele/Sviluppo/CApps/boost_1_52_0")
LINK_DIRECTORIES("G:/Daniele/Sviluppo/CApps/boost_1_52_0")

# Add compiler switches, add optimization later
add_definitions("-g -Wall -m64 -std=c++0x")
add_executable(${PROJECT_NAME} ${SRC_LIST})





