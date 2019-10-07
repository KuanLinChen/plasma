cmake_minimum_required(VERSION 2.6)


project(ultraMPP)

set (ultraMPP_VERSION_MAJOR 0)
set (ultraMPP_VERSION_MINOR 1)
set (ultraMPP_PATCH_VERSION 0)
set (ultraMPP_VERSION ${ultraMPP_VERSION_MAJOR}.${ultraMPP_VERSION_MINOR}.${ultraMPP_PATCH_VERSION})

set(CMAKE_DEBUG_POSTFIX -d)

list (APPEND CMAKE_MODULE_PATH "$ENV{ultraMPP_DIR}/cmake-modules")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -DUSEMPI")

# MPI is needed
find_package (MPI REQUIRED)

include_directories(${MPI_CXX_INCLUDE_PATH})

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated" )

find_package (PETSc REQUIRED)

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

include_directories( $ENV{ultraMPP_DIR}/include ${Boost_INCLUDE_DIR} ${PETSC_INCLUDES} $ENV{CGNS_INC_PATH} ${CMAKE_CURRENT_SOURCE_DIR})


link_directories( $ENV{ultraMPP_DIR}/lib $ENV{CGNS_LIB_PATH} $ENV{HDF5_LIB_PATH} )
	
add_subdirectory ( PlasmaFluid )

add_executable ( p3 main.cpp  )


target_link_libraries( p3 LINK_PUBLIC ${PETSC_LIBRARIES} cgns hdf5 libultraMPP.a libmetis.a ${MPI_LIBRARIES} 
${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${PETSC_LIBRARIES} PlasmaFluid )




