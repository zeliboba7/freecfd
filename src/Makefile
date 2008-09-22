CC=mpic++ 
#OPTS=-Wall -g
OPTS=-O3
EXE=fcfd

INCLUDES= -I$(CGNS_BASE)/include -I$(PARMETIS_BASE)/include 
LIBS= -L$(CGNS_BASE)/lib -L$(PARMETIS_BASE)/lib 

PETSC_LIBS = -L$(PETSC_DIR)/lib/$(PETSC_ARCH) -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc  -lnsl -lrt -L$(PETSC_DIR)/externalpackages/fblaslapack/$(PETSC_ARCH) -lflapack -lfblas -lgfortran -lX11 

PETSC_INCLUDE = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH)

SRC= \
	assemble_linear_system.cc \
	convective_face_flux.cc \
	diffusive_face_flux.cc \
	grid.cc grid_partition.cc grid_readers.cc  \
	initialize.cc initialize_linear_system.cc inputs.cc \
	main.cc mpi_functions.cc \
	petsc_functions.cc probes.cc \
	read_inputs.cc read_restart.cc \
	set_bcs.cc \
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
