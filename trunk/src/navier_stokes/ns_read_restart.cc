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
#include "ns.h"

void NavierStokes::read_restart(int restart_step,vector<vector<int> > &partitionMap) {

	string dirname="./restart/"+int2str(restart_step)+"/";
	string gs="."+int2str(gid+1);
	
	p.read_cell_data(dirname+"p"+gs,partitionMap);
	V.read_cell_data(dirname+"V"+gs,partitionMap);
	T.read_cell_data(dirname+"T"+gs,partitionMap);
	for (int c=0;c<grid[gid].cellCount;++c) rho.cell(c)=material.rho(p.cell(c),T.cell(c));
	update_boundaries();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();
	calc_limiter();

	return;
}


