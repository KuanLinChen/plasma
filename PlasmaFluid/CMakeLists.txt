cmake_minimum_required(VERSION 2.8)
project(PlasmaFluid)

set(CMAKE_DEBUG_POSTFIX -d)

# MPI is needed
find_package (MPI REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")



include_directories(${MPI_CXX_INCLUDE_PATH})
link_directories(${MPI_CXX_LIBRARIES})

include_directories( $ENV{ultraMPP_DIR}/include ${Boost_INCLUDE_DIR} ${PETSC_INCLUDES} $ENV{CGNS_INC_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

link_directories( $ENV{ultraMPP_DIR}/lib $ENV{CGNS_LIB_PATH} $ENV{HDF5_LIB_PATH} )


add_library ( 

PlasmaFluid 						PFM.hpp 
config_structure.cpp 		config_structure.hpp
domain_structure.cpp 		domain_structure.hpp
variable_structure.cpp 	variable_structure.hpp
chemistry-main.cpp			chemistry.h
solver_poisson.cpp 			solver_poisson.hpp
solver_drift_diffusion.cpp 	solver_drift_diffusion.hpp

post_structure.cpp 			post_structure.hpp
#variable_structure_NS.cpp 	variable_structure_NS.hpp
#solver_poisson_explicit.cpp solver_poisson_explicit.hpp
#solver_energy_density.cpp 	solver_energy_density.hpp
#solver_fluid_model.cpp 		solver_fluid_model.hpp
#solver_navier_stokes.cpp 	solver_navier_stokes.hpp

)

target_include_directories ( PlasmaFluid PUBLIC libultraMPP.a ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries ( PlasmaFluid LINK_PUBLIC ${MPI_LIBRARIES} )


