cmake_minimum_required(VERSION 2.6)


project(ultraMPP)

set (ultraMPP_VERSION_MAJOR 2)
set (ultraMPP_VERSION_MINOR 0)
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
find_package (Boost 1.5 COMPONENTS system filesystem REQUIRED)


#hdf5
set(temps $ENV{HDF5_LIB_PATH} )
if (temps EQUAL "")

else()
message ( "Adding hdf5 library")
message ("HDF5_library =  $ENV{HDF5_LIB_PATH}")
link_libraries( $ENV{HDF5_LIB_PATH} )

endif()




include_directories( $ENV{ultraMPP_DIR}/include ${Boost_INCLUDE_DIR} ${PETSC_INCLUDES} $ENV{CGNS_INC_PATH} ${CMAKE_CURRENT_SOURCE_DIR})


link_directories( $ENV{ultraMPP_DIR}/lib $ENV{CGNS_LIB_PATH} $ENV{HDF5_LIB_PATH} )
	


add_executable ( p3 interpolation.cpp  )


target_link_libraries( p3 LINK_PUBLIC ${PETSC_LIBRARIES} cgns hdf5 libultraMPP.a libmetis.a ${MPI_LIBRARIES} 
${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${PETSC_LIBRARIES} )

