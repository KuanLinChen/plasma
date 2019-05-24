cmake_minimum_required(VERSION 2.8)
project(system)

set(CMAKE_DEBUG_POSTFIX -d)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#add_definitions(-std=c++0x)

# MPI is needed
find_package (MPI REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
include_directories(${MPI_CXX_INCLUDE_PATH})
link_directories(${MPI_CXX_LIBRARIES})

# Boost
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
find_package (Boost 1.5 COMPONENTS system filesystem REQUIRED)
if (Boost_MINOR_VERSION LESS "50")
        set(BOOST_ROOT "/usr/local/")
        set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
        set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib64")
endif()

find_package (Boost 1.5 COMPONENTS system filesystem REQUIRED)

add_library (system STATIC sys_log.cpp sys_log.h uuid_functions.cpp uuid_functions.h tools.cpp tools.h )

target_include_directories (system PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(system LINK_PUBLIC ${MPI_LIBRARIES})

 
