CC=mpic++ 
#OPTS=-Wall -g
OPTS=-O3
EXE=fcfd

INCLUDES= -I$(CGNS_BASE)/include -I$(PARMETIS_BASE)/include 
LIBS= -L$(CGNS_BASE)/lib -L$(PARMETIS_BASE)/lib 

PETSC_LIBS = -L$(PETSC_DIR)/lib/$(PETSC_ARCH) -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc  -lnsl -lrt -L$(PETSC_DIR)/externalpackages/fblaslapack/$(PETSC_ARCH) -lflapack -lfblas -lgfortran -lX11 

PETSC_INCLUDE = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH)

SRC= \
	main.cc grid.cc vec3d.cc RoeFlux.cc  write_tec.cc write_vtk.cc writeTecplotMacro.cc \
        inputs.cc sparse.cc readInputs.cc initialize.cc setBCs.cc inviscid_flux.cc read_tec.cc \
        diffFlux.cc jac.cc

.SUFFIXES: .cc .o

OBJ=	$(SRC:.cc=.o)

.cc.o:
	$(CC) $(OPTS) $(INCLUDES) $(PETSC_INCLUDE) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(OPTS) -o $@ $(OBJ) $(LIBS) -lcgns -lparmetis -lmetis $(PETSC_LIBS)

clean:	
	rm -rf $(OBJ) $(EXE)
