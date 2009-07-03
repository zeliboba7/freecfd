/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

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
#include "rans.h"

// Updates the eddy viscosity stores at the cell centers
void RANS::update_cell_eddy_viscosity(void) {
	
	double arg1,arg2,arg3,F2;
	double mu=viscosity;
	double a1=0.31; // SST a1 value
	double turbulent_length_scale;
	for (int c=0;c<grid.cellCount;++c) {
		if (TURBULENCE_MODEL==SST) {
			arg1=2.*sqrt(cell[c].k+1.e-15)/(kepsilon.beta_star*cell[c].omega*grid.cell[c].closest_wall_distance);
			arg2=500.*mu/(grid.cell[c].rho*cell[c].omega
					*grid.cell[c].closest_wall_distance*grid.cell[c].closest_wall_distance);
			arg3=max(arg1,arg2);
			F2=tanh(arg3*arg3);
			cell[c].mu_t=a1*grid.cell[c].rho*cell[c].k/max(a1*cell[c].omega,cell[c].strainRate*F2);
			turbulent_length_scale=sqrt((cell[c].k+1.e-15))/max(cell[c].omega,cell[c].strainRate*F2/a1);
			turbulent_length_scale/=0.083;
		} else {
			cell[c].mu_t=grid.cell[c].rho*cell[c].k/cell[c].omega;
			turbulent_length_scale=sqrt(cell[c].k+1.e-15)/cell[c].omega;
			turbulent_length_scale/=0.083;
		}
		if (TURBULENCE_FILTER!=NONE) {
			double delta=turbulenceFilterSize;
			if (TURBULENCE_FILTER==LOCAL) delta*=grid.cell[c].lengthScale;
			cell[c].filterFunction=min(1.,delta/turbulent_length_scale);
			cell[c].mu_t*=cell[c].filterFunction;
		}
		
	}
	
	return;
	
} // end RANS::update_cell_eddy_viscosity

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



