#include <cmath>
#include<string>
#include "grid.h"
#include "bc.h"

extern Grid grid;
extern BC bc;

void diff_flux(double mu) {
	Vec3D tau_x, tau_y, tau_z;
	Vec3D gradUf, gradVf, gradWf,averageVel, areaVec;
	double diffFlux[4], factor;
	int parent, neighbor;
	for (int i=0; i<4; ++i) diffFlux[i]=0.;
	for (unsigned int f=0;f<grid.faceCount;++f) {
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
		areaVec=grid.face[f].normal*grid.face[f].area;
		averageVel=0.;
		for (unsigned int i=0;i<grid.face[f].cellContributions.indices.size();++i) {
			factor=grid.face[f].cellContributions.data[i];
			averageVel+=grid.cell[grid.face[f].cellContributions.indices[i]].v*factor;
		}
		if (grid.face[f].bc==-1) { //internal face
			gradUf=0.5* (grid.cell[parent].grad[1]+grid.cell[neighbor].grad[1]);
			gradVf=0.5* (grid.cell[parent].grad[2]+grid.cell[neighbor].grad[2]);
			gradWf=0.5* (grid.cell[parent].grad[3]+grid.cell[neighbor].grad[3]);
		} else {
			gradUf=grid.cell[parent].grad[1];
			gradVf=grid.cell[parent].grad[2];
			gradWf=grid.cell[parent].grad[3];
			if (bc.region[grid.face[f].bc].type=="slip") {
				averageVel-=averageVel.dot(grid.face[f].normal) *grid.face[f].normal;
			} else if (bc.region[grid.face[f].bc].type=="noslip") {
				averageVel=0.;
			} else if (bc.region[grid.face[f].bc].type=="inlet") {
				averageVel=bc.region[grid.face[f].bc].v;
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
		diffFlux[3]=tau_x.dot(averageVel) *areaVec.comp[0]+tau_y.dot(averageVel) *areaVec.comp[1]+tau_z.dot(averageVel) *areaVec.comp[2];
		for (int i = 1;i <=4;++i) {
			grid.cell[parent].flux[i] -= diffFlux[i-1];
			if (grid.face[f].bc==-1) {  // internal face
				grid.cell[neighbor].flux[i] += diffFlux[i-1];
			}
		}
	}

}
