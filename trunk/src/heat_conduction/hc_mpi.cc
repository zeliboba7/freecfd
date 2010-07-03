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
#include "hc.h"

void HeatConduction::mpi_init(void) {
	
	// Current processor number and the total number of processors
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	return;
}

void HeatConduction::mpi_update_ghost_primitives(void) {
	
	// Store the current time step values in the update array
	for (int g=0;g<grid[gid].ghostCount;++g) update.ghost(g)=T.ghost(g);
	
	// Send and receive ghost data
	T.mpi_update();

	// Get the difference between the old and new values
	for (int g=0;g<grid[gid].ghostCount;++g) update.ghost(g)=T.ghost(g)-update.ghost(g);

	return;
} 

void HeatConduction::mpi_update_ghost_gradients(void) {
	
	gradT.mpi_update();
	
	return;
} 

