project(mpFilterApp)
cmake_minimum_required(VERSION 2.6)
aux_source_directory(. SRC_LIST)
include_directories(".")

# SET MODULE PATH
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")

# INCLUDE BOOST AND MPI LIBS
find_package(Boost REQUIRED)
find_package(MPI REQUIRED)
find_package(METIS REQUIRED)
find_package(PARMETIS REQUIRED)

# INCLUDE LIBRARY PATH
include_directories(${Boost_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${METIS_INCLUDE_DIR})
include_directories("/usr/local/lib/")

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
add_definitions("-g -Wall -m64 -v -O0 -std=c++0x -fbounds-check -Wno-unused-variable")

add_executable(${PROJECT_NAME} ${SRC_LIST})
#
target_link_libraries(mpFilterApp ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES} parmetis)

