src = \
ns.cc \
ns_apply_bcs.cc \
ns_assemble_linear_system.cc \
ns_ausm_plus_up.cc \
ns_convective_face_flux.cc \
ns_diffusive_face_flux.cc \
ns_initialize_linear_system.cc \
ns_limiters.cc \
ns_mpi.cc \
ns_petsc_functions.cc \
ns_read_restart.cc \
ns_roe.cc \
ns_sd_slau.cc \
ns_set_bcs.cc \
ns_sources.cc \
ns_write_restart.cc 

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

ns : $(objects)
	ar rcs libns.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp ns.h ../../include
	cp ns_state_cache.h ../../include
	mv libns.a ../../lib
