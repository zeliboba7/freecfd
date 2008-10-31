/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

	This file is a part of Free CFD

	Free CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Free CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#ifndef MPI_FUNCTIONS_H
#define MPI_FUNCTIONS_H

#include <mpi.h>
#include <vector>
using namespace std;

#include "grid.h"

extern Grid grid;
extern int np, Rank;
extern IndexMaps maps;

	
void mpi_init(int argc, char *argv[]);
void mpi_handshake(void);
void mpi_update_ghost_primitives(void);
void mpi_update_ghost_gradients(void);
void mpi_update_ghost_limited_gradients(void);

struct mpiGhost {
	unsigned int globalId;
	double vars[8];
};
	
struct mpiGrad {
	unsigned int globalId;
	double grads[21];
};

#endif
