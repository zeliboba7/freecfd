CC=mpic++
 
ARCH_TYPE := $(shell uname -m | sed "s/ //g")
ARCH_FLAGS_i686 =  -march=pentium3 -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer
ARCH_FLAGS_ppc = -fast  -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer
ARCH_FLAGS_ia64 = -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer -mb-step
ARCH_FLAGS_x86_64 = -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-signaling-nans -fstrict-aliasing -fomit-frame-pointer

OPTS = -O3 $(ARCH_FLAGS_$(ARCH_TYPE))
#OPTS = -O3
#OPTS=-Wall -g

EXE=fcfd

INCLUDES= -I$(CGNS_BASE)/include -I$(PARMETIS_BASE)/include 
LIBS= -L$(CGNS_BASE)/lib -L$(PARMETIS_BASE)/lib 

PETSC_LIBS = -L$(PETSC_DIR)/lib/$(PETSC_ARCH) -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lnsl -lrt

PETSC_INCLUDE = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH)

SRC= \
	assemble_linear_system.cc \
	check_inputs.cc convective_face_flux.cc commons.cc \
	diffusive_face_flux.cc \
	eos.cc \
	flamelet.cc flamelet_table.cc flamelet_terms.cc\
	grid.cc grid_averaging_metrics.cc grid_create_elements.cc grid_partition.cc grid_readers.cc  \
	initialize.cc initialize_linear_system.cc inputs.cc \
	main.cc mpi_functions.cc \
	petsc_functions.cc probes.cc \
	read_inputs.cc read_restart.cc \
	set_bcs.cc sources.cc \
	rans.cc rans_terms.cc\
	update_eddy_viscosity.cc update_face_mdot.cc\
	write_output.cc write_restart.cc \
	vec3d.cc 

.SUFFIXES: .cc .o

OBJ=	$(SRC:.cc=.o)

.cc.o:
	$(CC) $(OPTS) $(INCLUDES) $(PETSC_INCLUDE) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(OPTS) -o $@ $(OBJ) $(LIBS) -lcgns -lparmetis -lmetis $(PETSC_LIBS)

clean:	
	rm -rf $(OBJ) $(EXE)
