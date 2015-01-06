PROJECT(mpFilterApp)
cmake_minimum_required(VERSION 2.6)
aux_source_directory(. SRC_LIST)
include_directories(".")

# SET MODULE PATH
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")

# INCLUDE BOOST AND MPI LIBS
find_package(Boost REQUIRED)
find_package(MPI REQUIRED)
find_package(METIS REQUIRED)
# find_package(PARMETIS)

# INCLUDE LIBRARY PATH
include_directories(${Boost_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${METIS_INCLUDE_DIR})
# include_directories("/usr/local/lib/")

# SET PLATFORM SPECIFIC OPTIONS
if(APPLE)
  set(CMAKE_C_COMPILER "/opt/local/bin/gcc")
  set(CMAKE_CXX_COMPILER "/opt/local/bin/g++")
  message(${MPI_C_FOUND})
  message(${MPI_CXX_FOUND})
  message(${MPI_C_COMPILER})
  message(${MPI_CXX_COMPILER})
  message(${MPI_INCLUDE_PATH})
  # INCLUDE DIRECTORIES FOR MPI
  # include_directories("/opt/local/include/mpich-gcc46")
  # link_directories("/opt/local/lib/mpich-gcc46")
endif()

# CHECK COMPILER
message( "USING COMPILER: " ${CMAKE_CXX_COMPILER_ID} )

# SET COMPILER FLAGS
add_definitions("-g -Wall -m64 -v -O0 -std=c++0x -Wno-unused-variable")

# CREATE EXECUTABLE
ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_LIST})

# LINK LIBRARIES
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${METIS_LIBRARIES})

