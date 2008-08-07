/************************************************************************
	
	Copyright 2007-2008 Emre Sozer & Patrick Clark Trizila

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
#include <cmath>
#include "grid.h"
#include "bc.h"
#include "petsc_functions.h"

extern Grid grid;
extern BC bc;
extern int rank;
extern double Gamma;

void roe_flux(const double qL[], const double qR[], double flux[]);

void viscous_jac(double mu) {

	int row,col;
	PetscScalar value;

	Vec3D tau_x, tau_y, tau_z;
	Vec3D gradUf, gradVf, gradWf, faceVel, areaVec;
	double viscousFlux[4], factor;
	int parent, neighbor,c;
	map<int,double>::iterator fit;
	map<int,Vec3D> stencil;
	map<int,Vec3D>::iterator git,sit;
				
	for (int i=0; i<4; ++i) viscousFlux[i]=0.;
	
	for (unsigned int f=0;f<grid.faceCount;++f) {
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
		areaVec=grid.face[f].normal*grid.face[f].area;

		// Use a map for the cells involved viscous fluxes and store contributions
		stencil.clear();
		
		for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
			c=(*fit).first;
				for (git=grid.cell[c].gradMap.begin();git!=grid.cell[c].gradMap.end(); git++ ) {
					if (stencil.find((*git).first)!=stencil.end()) {
						stencil[(*git).first]+=(*fit).second*(*git).second;
					} else {
						stencil.insert(pair<int,Vec3D>((*git).first,(*fit).second*(*git).second));
					}

						//cell[c].grad[0]+=(*git).second*cell[(*git).first].rho;
				} // end gradMap loop	
		}

		
		for (sit=stencil.begin();sit!=stencil.end(); sit++ ) {
			row=grid.cell[parent].globalId*5+1; // viscous_flux_x
			col=grid.cell[(*sit).first].globalId*5+1; //u
			// (viscous_flux_x due to u)=4/3*mu*du/dx*Ax+mu*du/dy*Ay+mu*du/dz*Az
			//                          =1/3*mu*du/dx*Ax+mu*(gradU.A)
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_x by u components
			value=1./3.*mu*(*sit).second.comp[0]*areaVec.comp[0]
					+mu*((*sit).second.dot(areaVec));
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================
			row=grid.cell[parent].globalId*5+1; // viscous_flux_x
			col=grid.cell[(*sit).first].globalId*5+2; //v
			// (viscous_flux_x due to v)=-2/3*mu*dv/dy*Ax+mu*dv/dx*Ay
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_x by v components
			value=-2./3.*mu*(*sit).second.comp[1]*areaVec.comp[0]
					+mu*(*sit).second.comp[0]*areaVec.comp[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================
			row=grid.cell[parent].globalId*5+1; // viscous_flux_x
			col=grid.cell[(*sit).first].globalId*5+3; //w
			// (viscous_flux_x due to w)=-2/3*mu*dw/dz*Ax+mu*dw/dx*Az
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_x by w components
			value=-2./3.*mu*(*sit).second.comp[2]*areaVec.comp[0]
					+mu*(*sit).second.comp[0]*areaVec.comp[2];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================
			row=grid.cell[parent].globalId*5+2; // viscous_flux_y
			col=grid.cell[(*sit).first].globalId*5+2; //v
			// (viscous_flux_y due to v)=mu*dv/dx*Ax+4/3*mu*dv/dy*Ay+mu*dv/dz*Az
			//                          =1/3*mu*dv/dy*Ay+mu*(gradV.A)
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_y by v components
			value=1./3.*mu*(*sit).second.comp[1]*areaVec.comp[1]
					+mu*((*sit).second.dot(areaVec));
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================
			row=grid.cell[parent].globalId*5+2; // viscous_flux_y
			col=grid.cell[(*sit).first].globalId*5+1; //u
			// (viscous_flux_y due to u)=-2/3*mu*du/dx*Ay+mu*du/dy*Ax
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_y by u components
			value=-2./3.*mu*(*sit).second.comp[0]*areaVec.comp[1]
					+mu*(*sit).second.comp[1]*areaVec.comp[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================
			row=grid.cell[parent].globalId*5+2; // viscous_flux_y
			col=grid.cell[(*sit).first].globalId*5+3; //w
			// (viscous_flux_y due to w)=-2/3*mu*dw/dz*Ay+mu*dw/dy*Az
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_y by w components
			value=-2./3.*mu*(*sit).second.comp[2]*areaVec.comp[1]
					+mu*(*sit).second.comp[1]*areaVec.comp[2];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================
			row=grid.cell[parent].globalId*5+3; // viscous_flux_z
			col=grid.cell[(*sit).first].globalId*5+3; //w
			// (viscous_flux_z due to w)=mu*dw/dx*Ax+mu*dw/dy*Ay+4/3*mu*dw/dz*Az
			//                          =1/3*mu*dw/dz*Az+mu*(gradV.A)
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_z by w components
			value=1./3.*mu*(*sit).second.comp[2]*areaVec.comp[2]
					+mu*((*sit).second.dot(areaVec));
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================
			row=grid.cell[parent].globalId*5+3; // viscous_flux_z
			col=grid.cell[(*sit).first].globalId*5+1; //u
			// (viscous_flux_z due to u)=-2/3*mu*du/dx*Az+mu*du/dz*Ax
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_z by u components
			value=-2./3.*mu*(*sit).second.comp[0]*areaVec.comp[2]
					+mu*(*sit).second.comp[2]*areaVec.comp[0];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================
			row=grid.cell[parent].globalId*5+3; // viscous_flux_z
			col=grid.cell[(*sit).first].globalId*5+2; //v
			// (viscous_flux_z due to v)=-2/3*mu*dv/dy*Az+mu*dv/dz*Ay
			// grad_contributions=(*sit).second
			// Hence total contribution to viscous_flux_z by v components
			value=-2./3.*mu*(*sit).second.comp[1]*areaVec.comp[2]
					+mu*(*sit).second.comp[2]*areaVec.comp[1];
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			//===================================================================

		} 
		
	}

	return;
} // end function


