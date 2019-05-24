cmake_minimum_required(VERSION 2.8)
project(metis)

set(GKLIB_PATH "${PROJECT_SOURCE_DIR}/GKlib" CACHE PATH "path to GKlib" FORCE )
set(SHARED FALSE CACHE BOOL "build a shared library" FORCE)

#message("GKLIB_PATH= ${GKLIB_PATH}")
#message("PROJECTPATH= ${PROJECT_SOURCE_DIR}")

if(MSVC)
  set(METIS_INSTALL FALSE)
else()
  set(METIS_INSTALL TRUE)
endif()

# Configure libmetis library.
if(SHARED)
  set(METIS_LIBRARY_TYPE SHARED)
else()
  set(METIS_LIBRARY_TYPE STATIC)
endif(SHARED)

include(${GKLIB_PATH}/GKlibSystem.cmake)
# Add include directories.
include_directories(${GKLIB_PATH})
include_directories(include)
# Recursively look for CMakeLists.txt in subdirs.
add_subdirectory("include")
add_subdirectory("libmetis")
#add_subdirectory("programs")
