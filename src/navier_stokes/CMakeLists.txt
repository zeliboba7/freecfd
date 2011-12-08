set (NAME ns)
set (SOURCES 
ns.cc
ns_convective_face_flux.cc     
ns_mpi.cc                      
ns_sd_slau.cc
ns_apply_bcs.cc                
ns_diffusive_face_flux.cc      
ns_petsc_functions.cc          
ns_set_bcs.cc
ns_assemble_linear_system.cc   
ns_read_restart.cc             
ns_sources.cc
ns_ausm_plus_up.cc             
ns_limiters.cc                 
ns_roe.cc                      
ns_stegger_warming.cc
ns_update_boundaries.cc
ns_vanleer.cc
ns_write_restart.cc
ns_time_terms.cc
)

add_library(${NAME} STATIC ${SOURCES} )

install (FILES ${NAME}.h DESTINATION include)
install (FILES lib${NAME}.a DESTINATION lib)
 
