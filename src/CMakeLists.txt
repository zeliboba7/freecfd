set (CMAKE_CXX_COMPILER icpc)
set (CMAKE_CXX_FLAGS "-O3")
#set (CMAKE_CXX_FLAGS "-fast")
#set (CGNS_INCLUDE_DIRS      /nobackup/esozer/codes/fcfd/deps/cgnslib_2.5/include)
#set (CGNS_LIBRARY_DIRS      /nobackup/esozer/codes/fcfd/deps/cgnslib_2.5/LINUX)
#set (PARMETIS_INCLUDE_DIRS  /nobackup/esozer/codes/fcfd/deps/ParMetis-3.1.1)
#set (PARMETIS_LIBRARY_DIRS  /nobackup/esozer/codes/fcfd/deps/ParMetis-3.1.1)
#set (PETSC_INCLUDE_DIRS     /nobackup/esozer/install/include)
#set (PETSC_LIBRARY_DIRS     /nobackup/esozer/install/lib)

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

set (CMAKE_INSTALL_PREFIX /home1/esozer/codes/freecfd)

# Intel(R) Compiler has its own library archiver,
# if you build libraries and do not use xiar,
# the Intel compiler will complain about invalid
# archives at the link phase.

# The Intel(R) archiver is "xiar" usually
# located in the same folder as the compiler,
# /opt/intel/Compiler/nn/n.nn/ia32/bin/xiar
# for example.

#FIND_PROGRAM(XIAR xiar)
#IF(XIAR)
#  SET(CMAKE_AR "${XIAR}")
#ENDIF(XIAR)
#MARK_AS_ADVANCED(XIAR)

# Intel(R) Compiler also comes with its own linker
# which provides a number of additional benefits when
# linking code compiled with the Intel(R) compiler.
# Again, usually in the same place as icc itself,

# /opt/intel/COmpiler/nn/n.nn/ia32/bin/xild per e.g

FIND_PROGRAM(XILD xild)
IF(XILD)
  SET(CMAKE_LINKER "${XILD}")
ENDIF(XILD)
MARK_AS_ADVANCED(XILD) 


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
set (EXTRA_LIBS parmetis metis cgns petsc mpi)

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
 
