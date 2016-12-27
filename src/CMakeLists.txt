cmake_minimum_required(VERSION 2.8.9)

#Project name
project(myMCS)

#Specify directory of header files 
include_directories(include)
#Select all the source files
file(GLOB SOURCES "*.cpp")
#Compile the library
add_library(${PROJECT_NAME} SHARED ${SOURCES})
#Compile the executable
add_executable(mcsWrap ${SOURCES})
#forcing to remove "lib" prefix
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
#forcing 32bit compilation
#set(CMAKE_C_FLAGS "-m32")
#set(CMAKE_CXX_FLAGS "-m32")

#forcing debug compilation
set(CMAKE_BUILD_TYPE DEBUG)


