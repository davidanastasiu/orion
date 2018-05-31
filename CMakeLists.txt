cmake_minimum_required(VERSION 2.8)
project(ORION)

if(MSVC)
  set(ORION_INSTALL FALSE)
else()
  set(ORION_INSTALL TRUE)
endif()

# Include directories
include_directories(src)

# Recursively look for CMakeLists.txt in subdirectories
add_subdirectory("src")

