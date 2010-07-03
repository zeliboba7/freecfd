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

void HeatConduction::set_bcs(void) {
	
	// Loop through each boundary condition region and apply sequentially
	int count=input.section("grid",gid).subsection("BC",0).count;
	
	for (int b=0;b<count;++b) {
		// Store the reference to current BC region
		Subsection &region=input.section("grid",gid).subsection("BC",b);
		
		if (region.get_double("T").is_found) {
			T.fixedonBC[b]=true; T.bcValue[b].resize(1);
			T.bc(b)=region.get_double("T");
			bc[gid][b].thermalType=FIXED_T;
		} else if (region.get_double("qdot").is_found) {
			bc[gid][b].thermalType=FIXED_Q;
			qdot.fixedonBC[b]=true; qdot.bcValue[b].resize(1);
			qdot.bc(b)=region.get_double("qdot");
		} else {
			bc[gid][b].thermalType=ADIABATIC;
		}
		
		// TODO: The following is done just to expand the arrays. 
		// If in the future, the bcValue arrays are by default full size, get rid of this
		// Loop through the interfaces
		int b;
		for (int g=0;g<grid.size();++g) {
			for (int i=0;i<interface[g].size();++i) {
				if (interface[g][i].donor_grid==gid) { // If this is a donor
					if(interface[g][i].donor_var=="qdot") qdot.bcValue[interface[g][i].donor_bc].resize(grid[gid].boundaryFaceCount[interface[g][i].donor_bc][grid[gid].Rank]);
					// Flowfield variables do not need to be allocated at the bc's
				} // if donor
				else if (interface[g][i].recv_grid==gid) { // If this is a receiver
					b=interface[g][i].recv_bc;
					if(interface[g][i].recv_var=="qdot") {
						qdot.bcValue[b].resize(grid[gid].boundaryFaceCount[b][grid[gid].Rank]);
						for (int j=1;j<qdot.bcValue[b].size();++j) qdot.bcValue[b][j]=qdot.bcValue[b][0];
					}
					if(interface[g][i].recv_var=="T") {
						T.bcValue[b].resize(grid[gid].boundaryFaceCount[b][grid[gid].Rank]);
						for (int j=1;j<T.bcValue[b].size();++j) T.bcValue[b][j]=T.bcValue[b][0];
					}
				} // if receiver
			}
		}
		

	}

	return;
}
