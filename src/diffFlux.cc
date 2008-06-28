#include <cmath>
#include<string>
#include "grid.h"
#include "bc.h"

extern Grid grid;
extern BC bc;

void diff_flux(double mu) {

	Vec3D tau_x, tau_y, tau_z;
	Vec3D gradUf, gradVf, gradWf, faceVel, areaVec;
	double diffFlux[4], factor;
	int parent, neighbor;
	map<int,double>::iterator fit;
	
	for (int i=0; i<4; ++i) diffFlux[i]=0.;
	
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
				faceVel=bc.region[grid.face[f].bc].v;
			}
		}
		
		tau_x.comp[0]=2./3.*mu* (2.*gradUf.comp[0]-gradVf.comp[1]-gradWf.comp[2]);
		tau_x.comp[1]=mu* (gradUf.comp[1]+gradVf.comp[0]);
		tau_x.comp[2]=mu* (gradUf.comp[2]+gradWf.comp[0]);
		tau_y.comp[0]=tau_x.comp[1];
		tau_y.comp[1]=2./3.*mu* (2.*gradVf.comp[1]-gradUf.comp[0]-gradWf.comp[2]);
		tau_y.comp[2]=mu* (gradVf.comp[2]+gradWf.comp[1]);
		tau_z.comp[0]=tau_x.comp[2];
		tau_z.comp[1]=tau_y.comp[2];
		tau_z.comp[2]=2./3.*mu* (2.*gradWf.comp[2]-gradUf.comp[0]-gradVf.comp[1]);

		diffFlux[0]=tau_x.dot(areaVec);
		diffFlux[1]=tau_y.dot(areaVec);
		diffFlux[2]=tau_z.dot(areaVec);
		diffFlux[3]=tau_x.dot(faceVel) *areaVec.comp[0]+tau_y.dot(faceVel) *areaVec.comp[1]+tau_z.dot(faceVel) *areaVec.comp[2];
		
		for (int i = 1;i <=4;++i) {
			grid.cell[parent].flux[i] -= diffFlux[i-1];
			if (grid.face[f].bc==-1) {  // internal face
				grid.cell[neighbor].flux[i] += diffFlux[i-1];
			}
		} // for i
		
	} // for face f

}
