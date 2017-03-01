cmake_minimum_required(VERSION 2.8.9)

#Project name
project(MCS)

#Specify directory of header files 
include_directories(include)
#Select all the source files
file(GLOB SOURCES "*.cpp")
#Compile the library
add_library(libMCS SHARED ${SOURCES})
#Compile the executable
add_executable(runMCS ${SOURCES})
#forcing to add "lib" prefix to library file
set_target_properties(libMCS PROPERTIES PREFIX "")
#forcing 32bit compilation
#set(CMAKE_C_FLAGS "-m32")
#set(CMAKE_CXX_FLAGS "-m32")

#forcing debug compilation
set(CMAKE_BUILD_TYPE DEBUG)


