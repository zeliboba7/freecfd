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

void HeatConduction::read_restart(int restart_step,vector<vector<int> > &partitionMap) {

	string dirname="./restart/"+int2str(restart_step)+"/";
	string gs="."+int2str(gid+1);
	
	T.read_cell_data(dirname+"T"+gs,partitionMap);

	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients(); 
}


