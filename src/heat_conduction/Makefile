src = \
hc.cc \
hc_assemble_linear_system.cc \
hc_diffusive_face_flux.cc \
hc_initialize_linear_system.cc \
hc_mpi.cc \
hc_petsc_functions.cc \
hc_read_restart.cc \
hc_set_bcs.cc \
hc_sources.cc \
hc_write_restart.cc 

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

ns : $(objects)
	ar rcs libhc.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp hc.h ../../include
	cp hc_state_cache.h ../../include
	mv libhc.a ../../lib
