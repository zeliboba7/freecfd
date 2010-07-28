/************************************************************************
 
 Copyright 2007-2010 Emre Sozer 
 
 Contact: emresozer@freecfd.com
 
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
#include "rans.h"

void RANS::mpi_init(void) {
	
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	return;
}

void RANS::mpi_update_ghost_primitives(void) {
	
	for (int g=0;g<grid[gid].ghostCount;++g) {
		update[0].ghost(g)=k.ghost(g);
		update[1].ghost(g)=omega.ghost(g);
	}
	
	k.mpi_update();
	omega.mpi_update();
	mu_t.mpi_update();

	for (int g=0;g<grid[gid].ghostCount;++g) {
		update[0].ghost(g)=k.ghost(g)-update[0].ghost(g);
		update[1].ghost(g)=omega.ghost(g)-update[1].ghost(g);
	}
	
	return;
} 

void RANS::mpi_update_ghost_gradients(void) {
	
	gradk.mpi_update();
	gradomega.mpi_update();

	return;
} 

