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

		row=grid.cell[parent].globalId*5+1;
		for (sit=stencil.begin();sit!=stencil.end(); sit++ ) {
			col=grid.cell[(*sit).first].globalId*5+1; //rho*u
			value=(*sit).second.comp[1]; // du/dy
			value*=mu*areaVec.comp[1]; //tau_xy.dot.areaVec
			MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
		} 
		
	}

// 		if (grid.face[f].bc>=0) {
// 			if (bc.region[grid.face[f].bc].type=="slip") {
// 				faceVel-=faceVel.dot(grid.face[f].normal)*grid.face[f].normal;
// 			} else if (bc.region[grid.face[f].bc].type=="noslip") {
// 				faceVel=0.;
// 			} else if (bc.region[grid.face[f].bc].type=="inlet") {
// 				faceVel=bc.region[grid.face[f].bc].v;
// 			}
// 		}
// 		
// 		tau_x.comp[0]=2./3.*mu* (2.*gradUf.comp[0]-gradVf.comp[1]-gradWf.comp[2]);
// 		tau_x.comp[1]=mu* (gradUf.comp[1]+gradVf.comp[0]);
// 		tau_x.comp[2]=mu* (gradUf.comp[2]+gradWf.comp[0]);
// 		tau_y.comp[0]=tau_x.comp[1];
// 		tau_y.comp[1]=2./3.*mu* (2.*gradVf.comp[1]-gradUf.comp[0]-gradWf.comp[2]);
// 		tau_y.comp[2]=mu* (gradVf.comp[2]+gradWf.comp[1]);
// 		tau_z.comp[0]=tau_x.comp[2];
// 		tau_z.comp[1]=tau_y.comp[2];
// 		tau_z.comp[2]=2./3.*mu* (2.*gradWf.comp[2]-gradUf.comp[0]-gradVf.comp[1]);
// 
// 		viscousFlux[0]=tau_x.dot(areaVec);
// 		viscousFlux[1]=tau_y.dot(areaVec);
// 		viscousFlux[2]=tau_z.dot(areaVec);
// 		viscousFlux[3]=tau_x.dot(faceVel)*areaVec.comp[0]+tau_y.dot(faceVel)*areaVec.comp[1]+tau_z.dot(faceVel) *areaVec.comp[2];
// 
// 		for (int i = 0;i <4;++i) {
// 			grid.cell[parent].flux[i+1] -= viscousFlux[i];
// 			if (grid.face[f].bc==-1) {  // internal face
// 				grid.cell[neighbor].flux[i+1] += viscousFlux[i];
// 			}
// 		} // for i
// 		
// 	} // for face f
// 	
// 	//MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);

	return;
} // end function


