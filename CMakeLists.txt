# cmake file to build MultiNest
# 
# Usage: ($ mkdir build)
#        $ cd build
#        $ cmake ..
#        ($ make clean)
#        $ make -j 4 
#
# After a successful build you can safely remove all files and directories from build/
#
# Notes:
#   - Prefer 'cmake ..' in the build/ directory above 'cmake' in the
#     parent directory, because the former does an out-of-source build.


project(MultiNest)

# cmake version shouldn't be too old

cmake_minimum_required(VERSION 2.8)


# Specify the location of the header files

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


# Specify the simulator source files to be compiled
              
file(GLOB sourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)


# Set the compiler flags
# -O3 = level 3 optimization; -Wall = enable all compiler's warning messages
# -fmessage-length=0 = each error message on a single line (no line-wrapping)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -fmessage-length=0 -pthread")


# Create an executable target

add_executable(multinest ${sourceFiles})


