# cmake file to build Diamonds
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
#
#   - You can change the default compiler by using (e.g. on Mac OSX)"
#     $ cmake -D CMAKE_CXX_COMPILER=/usr/bin/clang++ ..


project(Diamonds)

# cmake version shouldn't be too old

cmake_minimum_required(VERSION 2.8)


# Specify the location of the header files

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


# Specify the simulator source files to be compiled
              
file(GLOB sourceFiles ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)


# Set the compiler flags
# First those common to both gcc and clang:
# -O3   = level 3 optimization; 
# -Wall = enable all compiler's warning messages

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")

# Then those specific for the compiler, for using C++11 support.

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -std=c++11 -Wno-deprecated-register -fPIC")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fPIC")
endif()

# Create a shared library target

add_library(diamonds SHARED ${sourceFiles})
add_library(static-diamonds STATIC ${sourceFiles})

# Install the library in the lib/ folder


install(TARGETS diamonds LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)


file(MAKE_DIRECTORY pyDiamonds)
file(MAKE_DIRECTORY pyDiamonds/__init__.py "")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/python_wrapper/include)

add_subdirectory(pybind11)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBRARIES OFF)
SET(CMAKE_EXE_LINKER_FLAGS "-static")

add_library(pyDiamonds MODULE python_wrapper/PyDiamonds.cpp)
target_link_libraries(pyDiamonds PRIVATE pybind11::module)
target_link_libraries(pyDiamonds PRIVATE static-diamonds)
set_target_properties(pyDiamonds PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
        SUFFIX "${PYTHON_MODULE_EXTENSION}")

add_custom_command(TARGET pyDiamonds POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pyDiamonds> ${CMAKE_SOURCE_DIR}/pyDiamonds/pyDiamonds.so)


