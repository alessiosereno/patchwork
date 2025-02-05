# CMake project file for patchwork
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(patchwork LANGUAGES Fortran)

# Set the source directory
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Collect all Fortran source files in the src directory
file(GLOB SOURCES ${SRC_DIR}/*.f90)

# Create the executable
add_executable(patchwork ${SOURCES})

# Specify the Fortran standard (optional, adjust if needed)
set_target_properties(patchwork PROPERTIES Fortran_STANDARD 2008)

# Set directory for .mod files
set_target_properties(patchwork PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Set directory for executable output
set_target_properties(patchwork PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)