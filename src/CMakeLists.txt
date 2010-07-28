set (CMAKE_CXX_COMPILER mpic++)
set (CMAKE_CXX_FLAGS "-O3")
set (CGNS_INCLUDE_DIRS /Users/esozer/fcfd/deps/cgnslib_2.5/include)
set (CGNS_LIBRARY_DIRS /Users/esozer/fcfd/deps/cgnslib_2.5/lib)
set (PARMETIS_INCLUDE_DIRS /Users/esozer/fcfd/deps/ParMetis-3.1.1)
set (PARMETIS_LIBRARY_DIRS /Users/esozer/fcfd/deps/ParMetis-3.1.1)
set (PETSC_INCLUDE_DIRS /Users/esozer/fcfd/deps/petsc-3.1-p2/include /Users/esozer/fcfd/deps/petsc-3.1-p2/darwin9.8.0-c-debug/include)
set (PETSC_LIBRARY_DIRS /Users/esozer/fcfd/deps/petsc-3.1-p2/darwin9.8.0-c-debug/lib)

cmake_minimum_required(VERSION 2.6)

project(freecfd)
set (freecfd_VERSION_MAJOR 2)
set (freecfd_VERSION_MINOR 0)

list (APPEND CMAKE_MODULE_PATH "${freecfd_SOURCE_DIR}/CMake")

# Pass some CMake settings to source code through a header file
configure_file (
	"${PROJECT_SOURCE_DIR}/cmake_vars.h.in"
	"${PROJECT_BINARY_DIR}/cmake_vars.h"
)

set (CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../)

# add to the include search path
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/grid")
include_directories("${PROJECT_SOURCE_DIR}/heat_conduction")
include_directories("${PROJECT_SOURCE_DIR}/inputs")
include_directories("${PROJECT_SOURCE_DIR}/material")
include_directories("${PROJECT_SOURCE_DIR}/navier_stokes")
include_directories("${PROJECT_SOURCE_DIR}/polynomial")
include_directories("${PROJECT_SOURCE_DIR}/rans")
include_directories("${PROJECT_SOURCE_DIR}/utilities")
include_directories("${PROJECT_SOURCE_DIR}/variable")
include_directories("${PROJECT_SOURCE_DIR}/vec3d")

include_directories(${CGNS_INCLUDE_DIRS})
include_directories(${PARMETIS_INCLUDE_DIRS})
include_directories(${PETSC_INCLUDE_DIRS})
link_directories(${CGNS_LIBRARY_DIRS})
link_directories(${PARMETIS_LIBRARY_DIRS})
link_directories(${PETSC_LIBRARY_DIRS})

add_subdirectory(grid)
add_subdirectory(heat_conduction)
add_subdirectory(inputs)
add_subdirectory(material)
add_subdirectory(navier_stokes)
add_subdirectory(polynomial)
add_subdirectory(rans)
add_subdirectory(utilities)
add_subdirectory(variable)
add_subdirectory(vec3d)

set (FCFD_LIBS grid hc inputs material ns polynomial rans utilities variable vec3d)
set (EXTRA_LIBS cgns metis parmetis petsc f2cblas f2clapack)

#add the executable
set (SOURCES
bc_interface.cc      main.cc              read_restart.cc      time_step.cc         write_restart.cc
bc_interface_sync.cc read_inputs.cc       set_bcs.cc           write_output.cc
)
add_executable(freecfd ${SOURCES})

target_link_libraries (freecfd ${FCFD_LIBS} ${EXTRA_LIBS})

install (TARGETS freecfd RUNTIME DESTINATION bin)
 
