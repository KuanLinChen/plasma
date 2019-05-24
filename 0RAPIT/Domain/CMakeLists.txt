cmake_minimum_required(VERSION 2.8)
project(domain)

set(CMAKE_DEBUG_POSTFIX -d)

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

add_library (domain domain.cpp domain.h element.cpp element.h element_indexing.cpp cgns.cpp cgns.h fieldview.cpp fieldview.h gmsh.cpp gmsh.h preprocessingtool.cpp preprocessingtool.h preprocessing.h preprocessing.cpp calculate_dA_dL.cpp calculate_face_sign.cpp calculate_face_coefficient.cpp calculate_cell_face_area.cpp facetagging.cpp surface.cpp surface.h cellbycell_tracking.cpp cellbycell_tracking.h connection.h connection.cpp structured_mesh.cpp pti.cpp pti.h gradient.cpp )

#include_directories (${metis_SOURCE_DIR}/include ${system_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include)

target_include_directories (domain PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR} )

target_link_libraries(domain LINK_PUBLIC ${MPI_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} cgns hdf5 metis system)

#add_executable (o_domain_test o_domain_test.cpp domain.h)
#target_link_libraries(o_domain_test LINK_PUBLIC domain)
