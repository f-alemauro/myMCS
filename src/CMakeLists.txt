cmake_minimum_required(VERSION 2.8.9)

#Project name
project(MCS)

#Specify directory of header files 
include_directories(include)
#Select all the source files
file(GLOB SOURCES "*.cpp")

#Compile the library
if (WIN32)
	SET(LIBNAME "libMCS")
elseif (UNIX AND NOT APPLE)
	if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
		SET(LIBNAME "libMCS-x64")
	else( CMAKE_SIZEOF_VOID_P EQUAL 4 )
		SET(LIBNAME "libMCS-x32")
	endif ()
endif ()


add_library(${LIBNAME} SHARED ${SOURCES})
#Compile the executable
add_executable(runMCS ${SOURCES})
#forcing to add "lib" prefix to library file
set_target_properties(${LIBNAME} PROPERTIES PREFIX "")
#forcing 32bit compilation
#set(CMAKE_C_FLAGS "-m32")
#set(CMAKE_CXX_FLAGS "-m32")

#forcing debug compilation
set(CMAKE_BUILD_TYPE DEBUG)


