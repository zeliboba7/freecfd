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
#include<string>
#include "grid.h"
#include "bc.h"
#include "petsc_functions.h"

extern Grid grid;
extern BC bc;

void viscous_flux(double mu) {

	Vec3D tau_x, tau_y, tau_z;
	Vec3D gradUf, gradVf, gradWf, faceVel, areaVec;
	double viscousFlux[4];
	int parent, neighbor;
	map<int,double>::iterator fit;
	Vec3D p2n; // parent to neighbor distance
	Vec3D p2nUnit;
	double p2nNorm;
	Vec3D neighborV;

	PetscInt row;
	PetscScalar value;
	
	for (int i=0; i<4; ++i) viscousFlux[i]=0.;
	
	for (unsigned int f=0;f<grid.faceCount;++f) {
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
		areaVec=grid.face[f].normal*grid.face[f].area;

		// Find face averaged variables
		faceVel=0.; gradUf=0.; gradVf=0.; gradWf=0.;
		for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
			if ((*fit).first>=0) { // if contribution is coming from a real cell
				faceVel+=(*fit).second*grid.cell[(*fit).first].v;
				gradUf+=(*fit).second*grid.cell[(*fit).first].grad[1];
				gradVf+=(*fit).second*grid.cell[(*fit).first].grad[2];
				gradWf+=(*fit).second*grid.cell[(*fit).first].grad[3];
			} else { // if contribution is coming from a ghost cell
				faceVel+=(*fit).second*grid.ghost[-1*((*fit).first+1)].v;
				gradUf+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[1];
				gradVf+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[2];
				gradWf+=(*fit).second*grid.ghost[-1*((*fit).first+1)].grad[3];
			}
		}


		
		if (grid.face[f].bc>=0) {
			if (bc.region[grid.face[f].bc].type=="slip") {
				faceVel-=faceVel.dot(grid.face[f].normal)*grid.face[f].normal;
			} else if (bc.region[grid.face[f].bc].type=="noslip") {
				faceVel=0.;
			} else if (bc.region[grid.face[f].bc].type=="inlet") {
				//cout << parent << "\t" << grid.cell[parent].v << "\t" << grid.cell[parent].grad[1] << endl;
				faceVel=bc.region[grid.face[f].bc].v;
			}
			p2n=(grid.face[f].centroid-grid.cell[parent].centroid).dot(grid.face[f].normal);
			neighborV=faceVel;
			
		} else {
			p2n=grid.cell[neighbor].centroid-grid.cell[parent].centroid;
			neighborV=grid.cell[neighbor].v;
		}

		p2nNorm=fabs(p2n);
		p2nUnit=p2n/p2nNorm;
		// Replace the face normal directional derivative
// 		gradUf-=gradUf.dot(p2nUnit)*p2nUnit;
// 		gradUf+=((neighborV.comp[0]-grid.cell[parent].v.comp[0])/p2nNorm)*p2nUnit;
// 
// 		gradVf-=gradVf.dot(p2nUnit)*p2nUnit;
// 		gradVf+=((neighborV.comp[1]-grid.cell[parent].v.comp[1])/p2nNorm)*p2nUnit;
// 
// 		gradWf-=gradWf.dot(p2nUnit)*p2nUnit;
// 		gradWf+=((neighborV.comp[2]-grid.cell[parent].v.comp[2])/p2nNorm)*p2nUnit;

		
		tau_x.comp[0]=2./3.*mu* (2.*gradUf.comp[0]-gradVf.comp[1]-gradWf.comp[2]);
		tau_x.comp[1]=mu* (gradUf.comp[1]+gradVf.comp[0]);
		tau_x.comp[2]=mu* (gradUf.comp[2]+gradWf.comp[0]);
		tau_y.comp[0]=tau_x.comp[1];
		tau_y.comp[1]=2./3.*mu* (2.*gradVf.comp[1]-gradUf.comp[0]-gradWf.comp[2]);
		tau_y.comp[2]=mu* (gradVf.comp[2]+gradWf.comp[1]);
		tau_z.comp[0]=tau_x.comp[2];
		tau_z.comp[1]=tau_y.comp[2];
		tau_z.comp[2]=2./3.*mu* (2.*gradWf.comp[2]-gradUf.comp[0]-gradVf.comp[1]);

		viscousFlux[0]=tau_x.dot(areaVec);
		viscousFlux[1]=tau_y.dot(areaVec);
		viscousFlux[2]=tau_z.dot(areaVec);
		viscousFlux[3]=tau_x.dot(faceVel)*areaVec.comp[0]+tau_y.dot(faceVel)*areaVec.comp[1]+tau_z.dot(faceVel) *areaVec.comp[2];

		for (int i = 0;i <4;++i) {
			grid.cell[parent].flux[i+1] -= viscousFlux[i];
			if (grid.face[f].bc==-1) {  // internal face
				grid.cell[neighbor].flux[i+1] += viscousFlux[i];
			}
		} // for i
		
// 		for (int i = 0;i <4;++i) {
// 			grid.cell[parent].flux[i+1] -= viscousFlux[i];
// 			if (grid.face[f].bc==-1) {  // internal face
// 				grid.cell[neighbor].flux[i+1] += viscousFlux[i];
// 			}
// 		} // for i
		
		for (int i = 0;i <4;++i) {
			row=grid.cell[parent].globalId*5+i+1;
			value=viscousFlux[i];
			VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			if (grid.face[f].bc==-1) { // TODO what if a ghost face??
				row=grid.cell[neighbor].globalId*5+i+1;
				value*=-1.;
				VecSetValues(rhs,1,&row,&value,ADD_VALUES);
			}
		}
		
	} // for face f

}
