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
#include "turbulence.h"

// Updates the eddy viscosity stores at the faces
void Turbulence::update_eddy_viscosity(void) {
	
	map<int,double>::iterator fit;
	double faceK,faceOmega,faceRho,mu_t;
	unsigned int parent;
	
	for (unsigned int f=0;f<grid.faceCount;++f) {
		parent=grid.face[f].parent;
		// Find face averaged variables
		faceRho=0.; faceK=0.; faceOmega=0.;
		for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
			if ((*fit).first>=0) { // if contribution is coming from a real cell
				faceK+=(*fit).second*cell[(*fit).first].k;
				faceOmega+=(*fit).second*cell[(*fit).first].omega;
				faceRho+=(*fit).second*grid.cell[(*fit).first].rho;
			} else { // if contribution is coming from a ghost cell
				faceK+=(*fit).second*ghost[-1*((*fit).first+1)].k;
				faceOmega+=(*fit).second*ghost[-1*((*fit).first+1)].omega;
				faceRho+=(*fit).second*grid.ghost[-1*((*fit).first+1)].rho;
			}
		}
		
		faceK=max(faceK,kLowLimit);
		faceOmega=max(faceOmega,omegaLowLimit);
		
		if (grid.face[f].bc>=0) { // boundary face
			if (bc.region[grid.face[f].bc].type==NOSLIP) {
				mu_t=0.;
			} else if (bc.region[grid.face[f].bc].type==SYMMETRY) {
				mu_t=grid.cell[parent].rho*cell[parent].k/cell[parent].omega;
			} else if (bc.region[grid.face[f].bc].type==SLIP) {
				mu_t=faceRho*faceK/faceOmega;
			} else if (bc.region[grid.face[f].bc].type==INLET) {
				mu_t=faceRho*bc.region[grid.face[f].bc].k/bc.region[grid.face[f].bc].omega;
			} else if (bc.region[grid.face[f].bc].type==OUTLET) {
				mu_t=faceRho*faceK/faceOmega;
			}
		}
		
		face[f].mu_t=min(fabs(mu_t),viscosityRatioLimit*viscosity);
		
	}
	
} // end Turbulence::update_eddy_viscosity




