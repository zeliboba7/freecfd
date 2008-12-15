#########################################################################
#	
#	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila
#
#	Contact: emresozer@freecfd.com , ptrizila@freecfd.com
#
#	This file is a part of Free CFD
#
#	Free CFD is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   any later version.
#
#   Free CFD is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   For a copy of the GNU General Public License,
#   see <http://www.gnu.org/licenses/>.
#
#########################################################################
CC=mpic++ 
#OPTS=-Wall -g
OPTS=-O3
EXE=fcfd

INCLUDES= -I$(CGNS_BASE)/include -I$(PARMETIS_BASE)/include 
LIBS= -L$(CGNS_BASE)/lib -L$(PARMETIS_BASE)/lib 

PETSC_LIBS = -L$(PETSC_DIR)/lib/$(PETSC_ARCH) -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lnsl -lrt

PETSC_INCLUDE = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/bmake/$(PETSC_ARCH)

SRC= \
	assemble_linear_system.cc \
	commons.cc convective_face_flux.cc check_inputs.cc\
	diffusive_face_flux.cc \
	grid.cc grid_create_elements.cc grid_partition.cc grid_readers.cc  \
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
