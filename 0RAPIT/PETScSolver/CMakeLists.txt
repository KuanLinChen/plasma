cmake_minimum_required(VERSION 2.8)
project(PETScSolver)

set(CMAKE_DEBUG_POSTFIX -d)

find_package (PETSc REQUIRED)
#include_directories (${PETSC_INCLUDES} ${CMAKE_SOURCE_DIR}/include)
add_definitions (${PETSC_DEFINITIONS})

add_library (PETScSolver PETSc_solver.cpp PETSc_solver.h)

target_include_directories (PETScSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${PETSC_INCLUDES} )

target_link_libraries ( PETScSolver LINK_PUBLIC system ${PETSC_LIBRARIES} ${MPI_LIBRARIES}) 