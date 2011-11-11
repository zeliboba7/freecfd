set (CMAKE_CXX_COMPILER mpic++)
set (CMAKE_CXX_FLAGS "-O3")
set (CGNS_INCLUDE_DIRS      /nobackup/esozer/codes/fcfd/deps/cgnslib_2.5/include)
set (CGNS_LIBRARY_DIRS      /nobackup/esozer/codes/fcfd/deps/cgnslib_2.5/LINUX)
set (PARMETIS_INCLUDE_DIRS  /nobackup/esozer/codes/fcfd/deps/ParMetis-3.1.1)
set (PARMETIS_LIBRARY_DIRS  /nobackup/esozer/codes/fcfd/deps/ParMetis-3.1.1)
set (PETSC_INCLUDE_DIRS     /nobackup/esozer/install/include)
set (PETSC_LIBRARY_DIRS     /nobackup/esozer/install/lib)

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

set (CMAKE_INSTALL_PREFIX /nobackup/esozer/install)

# add to the include search path
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/grid")
include_directories("${PROJECT_SOURCE_DIR}/heat_conduction")
include_directories("${PROJECT_SOURCE_DIR}/inputs")
include_directories("${PROJECT_SOURCE_DIR}/interpolate")
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
add_subdirectory(interpolate)
add_subdirectory(material)
add_subdirectory(navier_stokes)
add_subdirectory(polynomial)
add_subdirectory(rans)
add_subdirectory(utilities)
add_subdirectory(variable)
add_subdirectory(vec3d)

set (FCFD_LIBS grid hc inputs interpolate material ns polynomial rans utilities variable vec3d)
set (EXTRA_LIBS parmetis metis cgns petsc)

#add the executable
set (SOURCES
bc_interface.cc
curvilinear_grad_map.cc
face_interpolation_weights.cc
gradient_maps.cc
lsqr_grad_map.cc
main.cc
node_interpolation_weights.cc
read_restart.cc
time_step.cc
write_restart.cc			
bc_interface_sync.cc	
read_inputs.cc       
set_bcs.cc 
write_surface_output.cc	
write_volume_output.cc
)

add_executable(freecfd ${SOURCES})

target_link_libraries (freecfd ${FCFD_LIBS} ${EXTRA_LIBS})

install (TARGETS freecfd RUNTIME DESTINATION bin)
 
