cmake_minimum_required(VERSION 2.8)
project(ORION)

set(GKLIB_PATH "GKlib" CACHE PATH "path to GKlib")

# Add GKlib's cmake settings
include(${GKLIB_PATH}/GKlibSystem.cmake)

# Source files
file(GLOB orion_sources *.c)

# Include directories
include_directories(${GKLIB_PATH})
include_directories(.)

# Build
add_executable(orion ${GKlib_sources} ${orion_sources})

foreach(prog orion)
  target_link_libraries(${prog} m)
endforeach(prog)

# Install
if(ORION_INSTALL)
  install(TARGETS orion 
    RUNTIME DESTINATION bin)
endif()

# Try to find subversion revision.
set(SVNREV "")
file(TO_NATIVE_PATH ${PROJECT_SOURCE_DIR}/../../ svn_dir)
message("  -- SVN dir set to ${svn_dir}")
if(IS_DIRECTORY "${svn_dir}.svn")
  include(FindSubversion)
  if(Subversion_FOUND)
    Subversion_WC_INFO(${svn_dir} orion)
    set(SVNREV ${orion_WC_REVISION})
  endif(Subversion_FOUND)
endif()
add_definitions(-DSVNINFO="${SVNREV}")
