set (NAME hc)
set (SOURCES 
hc.cc                          hc_diffusive_face_flux.cc      hc_mpi.cc                      hc_read_restart.cc             hc_sources.cc
hc_assemble_linear_system.cc   hc_initialize_linear_system.cc hc_petsc_functions.cc          hc_set_bcs.cc                  hc_write_restart.cc
)

add_library(${NAME} STATIC ${SOURCES} )

install (FILES ${NAME}.h DESTINATION include)
install (FILES lib${NAME}.a DESTINATION lib)
 
