cmake_minimum_required(VERSION 2.8)
project(chemistry)

set(CMAKE_DEBUG_POSTFIX -d)

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

add_library (chemistry chemistry.cpp chemistry.h reaction.cpp reaction.h species.cpp species.h)

target_include_directories (chemistry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/include)
target_link_libraries(chemistry LINK_PUBLIC config++)

# LibConfig
# find_package (Config++ REQUIRED)
