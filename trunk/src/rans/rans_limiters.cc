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

void RANS::calc_limiter(void) {
	if(order==FIRST) {
		// Do nothing
	} else {
		if (limiter_function==NONE) {
			// Do nothing
		} else {
			// Not working for now
			//barth_jespersen_limiter();
		}
	} 
	return;
}

void RANS::barth_jespersen_limiter(void) {
	
	int neighbor,g;
	double phi[2];
	double delta;
	double max_val[2], min_val[2];
	
	for (int c=0;c<grid[gid].cellCount;++c) {
		
		for (int i=0;i<2;++i) phi[i]=1.;
		
		max_val[0]=min_val[0]=k.cell(c);
		max_val[1]=min_val[1]=omega.cell(c);
		
		// First loop through face neighbors to find max values
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a bouddary face
				c==grid[gid].cellFace(c,cf).parent ? neighbor=grid[gid].cellFace(c,cf).neighbor : neighbor=grid[gid].cellFace(c,cf).parent;
				if (neighbor>=0) { // real cell
					max_val[0]=max(max_val[0],k.cell(neighbor));
					max_val[1]=max(max_val[1],omega.cell(neighbor));
					
					min_val[0]=min(min_val[0],k.cell(neighbor));
					min_val[1]=min(min_val[1],omega.cell(neighbor));
				} else { // neighbor cell is a ghost (partition interface)
					neighbor=-1*neighbor-1;
					max_val[0]=max(max_val[0],k.ghost(neighbor));
					max_val[1]=max(max_val[1],omega.ghost(neighbor));
					
					min_val[0]=min(min_val[0],k.ghost(neighbor));
					min_val[1]=min(min_val[1],omega.ghost(neighbor));
				}
			}
		} // end face loop
				
		// Second loop through face neigbors to calculate min limiter
		for (int cf=0;cf<grid[gid].cell[c].faceCount;++cf) {
			if (grid[gid].cellFace(c,cf).bc<0) { // If not a boundary face
				for (int var=0;var<2;++var) {
					if (var==0) {
						delta=gradk.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
						if (delta>0.) phi[var]=min(phi[var],min(1.,(max_val[var]-k.cell(c))/delta));
						else if (delta<0.) phi[var]=min(phi[var],min(1.,(min_val[var]-k.cell(c))/delta));
					} else if (var==1) {
						delta=gradomega.cell(c).dot(grid[gid].cellFace(c,cf).centroid-grid[gid].cell[c].centroid);
						if (delta>0.) phi[var]=min(phi[var],min(1.,(max_val[var]-omega.cell(c))/delta));
						else if (delta<0.) phi[var]=min(phi[var],min(1.,(min_val[var]-omega.cell(c))/delta));
					}
				}
			}
		}

		for (int var=0;var<2;++var) limiter[var].cell(c)=phi[var];
		
	} // end cell loop
	
	for (int var=0;var<2;++var) limiter[var].mpi_update();	
	
	return;
}
