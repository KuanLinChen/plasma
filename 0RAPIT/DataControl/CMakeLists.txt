cmake_minimum_required(VERSION 2.8)
project(datacontrol)

set(CMAKE_DEBUG_POSTFIX -d)

add_library (datacontrol datamanager.cpp datamanager.h dataexchanger.cpp dataexchanger_functions.cpp dataexchanger.h datapool.cpp datapool.h scalar.cpp scalar.h)

target_include_directories (datacontrol PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include )

target_link_libraries(datacontrol LINK_PUBLIC ${MPI_LIBRARIES} domain)

