set (NAME rans)
set (SOURCES 
rans.cc                       rans_petsc_functions.cc       rans_set_bcs.cc               rans_update_eddy_viscosity.cc
rans_mpi.cc                   rans_read_restart.cc          rans_terms.cc                 rans_write_restart.cc
)

add_library(${NAME} STATIC ${SOURCES} )

install (FILES ${NAME}.h DESTINATION include)
install (FILES lib${NAME}.a DESTINATION lib)
 
