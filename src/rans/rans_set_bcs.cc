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

void RANS::set_bcs(void) {
	
	// Loop through each boundary condition region and apply sequentially
	int count=input.section("grid",gid).subsection("BC",0).count;
	
	for (int b=0;b<count;++b) {
		// Store the reference to current BC region
		Subsection &region=input.section("grid",gid).subsection("BC",b);
		string type=region.get_string("type");
		string kind=region.get_string("kind");
		double intensity,viscRatio;
		// Assign specified values
		intensity=region.get_double("turbulenceintensity");
		viscRatio=region.get_double("eddyviscosityratio");

		if (type=="inlet") {
			if (!region.get_double("turbulenceintensity").is_found) {
				cerr << "[E] Turbulence intensity needs to be specified in grid=" << gid+1 << " boundary condition BC_" << b+1 << endl;
				exit(1);
			}
			if (!region.get_double("eddyviscosityratio").is_found) {
				cerr << "[E] Eddy viscosity ratio needs to be specified in grid=" << gid+1 << " boundary condition BC_" << b+1 << endl;
				exit(1);
			}
			k.fixedonBC[b]=true; omega.fixedonBC[b]=true;
			k.bcValue[b].resize(1); omega.bcValue[b].resize(1);

			// TODO: This could be done better
			// Find a face index that is on this boundary (if any on this partition)
			int fid=0;
			for (int f=0;f<grid[gid].faceCount;++f) {
				if (grid[gid].face[f].bc==b) {
					fid=f;
					break;
				}
			}
			double rho,T;
			rho=ns[gid].rho.face(fid);
			T=ns[gid].T.face(fid);
			k.bc(b)=intensity*fabs(ns[gid].V.face(fid));
			k.bc(b)*=1.5*k.bc(b);
			k.bc(b)=max(kLowLimit,k.bc(b));
	 	  	omega.bc(b)=rho*k.bc(b)/(viscRatio*material.viscosity(T));
			mu_t.fixedonBC[b]=true; mu_t.bcValue[b].resize(1);
			mu_t.bc(b)=viscRatio*material.viscosity(T);
		} else if (type=="wall") {
			if (kind!="slip") {
				mu_t.fixedonBC[b]=true; mu_t.bcValue[b].resize(1);
				mu_t.bc(b)=0.;		
			}
		}
		
		 yplus.bcValue[b].resize(grid[gid].boundaryFaceCount[b][grid[gid].Rank]);
		
	}
		
	return;
}

