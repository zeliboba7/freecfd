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

void RANS::update_eddy_viscosity(void) {
	
	double arg1,arg2,arg3,F2;
	double mu;
	double a1=0.31; // SST a1 value
	double turbulent_length_scale;

	for (int c=0;c<grid[gid].cellCount;++c) {
		mu=ns[gid].material.viscosity(ns[gid].T.cell(c));
		if (model==SST) {
			arg1=2.*sqrt(k.cell(c)+1.e-15)/(kepsilon.beta_star*omega.cell(c)*grid[gid].cell[c].closest_wall_distance);
			arg2=500.*mu/(ns[gid].rho.cell(c)*omega.cell(c)
					*grid[gid].cell[c].closest_wall_distance*grid[gid].cell[c].closest_wall_distance);
			arg3=max(arg1,arg2);
			F2=tanh(arg3*arg3);
			mu_t.cell(c)=a1*ns[gid].rho.cell(c)*k.cell(c)/max(a1*omega.cell(c),strainRate.cell(c)*F2);
			turbulent_length_scale=sqrt((k.cell(c)+1.e-15))/max(omega.cell(c),strainRate.cell(c)*F2/a1);
			turbulent_length_scale/=0.083;
		} else {
			mu_t.cell(c)=ns[gid].rho.cell(c)*k.cell(c)/omega.cell(c);
			turbulent_length_scale=sqrt(k.cell(c)+1.e-15)/omega.cell(c);
			turbulent_length_scale/=0.083;
		}
	}
	
	return;
	
} 

/*
// Updates the eddy viscosity stores at the faces
void RANS::update_face_eddy_viscosity(void) {
	
	// Update face center values
	map<int,double>::iterator fit;
	int parent;
	
	for (int f=0;f<grid.faceCount;++f) {
		parent=grid.face[f].parent;
		// Find face averaged value
		face[f].mu_t=0.;
		for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
			if ((*fit).first>=0) { // if contribution is coming from a real cell
				face[f].mu_t+=(*fit).second*cell[(*fit).first].mu_t;
			} else { // if contribution is coming from a ghost cell
				face[f].mu_t+=(*fit).second*ghost[-1*((*fit).first+1)].mu_t;
			}
		}
		
		face[f].mu_t=max(0.,face[f].mu_t);
		
		// Correct for boundaries
		if (grid.face[f].bc>=0) { // boundary face
			if (bc.region[grid.face[f].bc].type==NOSLIP) {
				face[f].mu_t=0.;
			} else if (bc.region[grid.face[f].bc].type==SYMMETRY) {
				face[f].mu_t=cell[parent].mu_t;
			}
		}
		
	}
	
	return;
	
} // end RANS::update_face_eddy_viscosity

*/


