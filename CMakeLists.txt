project(mpFilterApp)
cmake_minimum_required(VERSION 2.6)
aux_source_directory(. SRC_LIST)
include_directories(".")
# INCLUDE BOOST LIB
find_package(Boost 1.30 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
# CHECK COMPILER
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  message( using Clang )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  message( using GCC )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  message( using Intel C++ )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message( using Visual Studio C++ )
endif()
#
#add_definitions("-g -Wall -m64 -O0 -std=c++0x")
add_definitions("-g -Wall -m64 -O0")
add_executable(${PROJECT_NAME} ${SRC_LIST})
#
target_link_libraries(mpFilterApp ${Boost_LIBRARIES})
