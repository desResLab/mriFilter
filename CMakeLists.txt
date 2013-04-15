project(mpFilterApp)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
include_directories(".")

# Windows Compilation
#message("Including Boost Folders...")
#INCLUDE_DIRECTORIES("G:/Daniele/Sviluppo/CApps/boost_1_52_0")
#LINK_DIRECTORIES("G:/Daniele/Sviluppo/CApps/boost_1_52_0")

add_definitions("-g -Wall -m64 -O0 -std=c++0x")
add_executable(${PROJECT_NAME} ${SRC_LIST})
