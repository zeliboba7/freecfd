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

void inviscid_jac(void) {

	double rhoL,rhoR,pL,pR;
	Vec3D faceTangent1, faceTangent2, qVL, qVR, deltaV;
	double qNL[5], qNR[5], fluxNormal[5], flux[5], fluxPlus[5], deltaL[5], deltaR[5];
	double Mach;
	
	unsigned int parent,neighbor,f;
	double temp[5];
	int row,col;
	PetscScalar value;

	deltaL[0]=deltaR[0]=sqrt(std::numeric_limits<double>::epsilon());
	for (int i=1;i<5;++i) {
		deltaL[i]=deltaR[i]=deltaL[i-1];
	}
	deltaV=0.;

	for (unsigned int f=0;f<grid.faceCount;++f) {
		parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;

		// Take 2 nodes of the face to find a vector in the face plane
		faceTangent1= (grid.face[f].node(1)-grid.face[f].node(0));
		faceTangent1/=fabs(faceTangent1);
		// Cross the tangent vector with the normal vector to get the second tangent
		faceTangent2= (grid.face[f].normal).cross(faceTangent1);

		rhoL=grid.cell[parent].rho;
		pL=grid.cell[parent].p;
		qVL=rhoL*grid.cell[parent].v;

		qNL[0]=rhoL;
		qNL[1]=qVL.dot(grid.face[f].normal);
		qNL[2]=qVL.dot(faceTangent1);
		qNL[3]=qVL.dot(faceTangent2);
		qNL[4]=0.5*(qVL.dot(qVL))/rhoL +pL/ (Gamma - 1.);

		if (grid.face[f].bc==-1) { // means a real internal face

			rhoR=grid.cell[neighbor].rho;
			pR=grid.cell[neighbor].p;
			qVR=rhoR*grid.cell[neighbor].v;

			qNR[0]=rhoR;
			qNR[1]=qVR.dot(grid.face[f].normal);
			qNR[2]=qVR.dot(faceTangent1);
			qNR[3]=qVR.dot(faceTangent2);
			qNR[4]=0.5*(qVR.dot(qVR))/rhoR +pR/ (Gamma - 1.);

		} else if (grid.face[f].bc>=0) {
			// Find face averaged variables
			map<int,double>::iterator fit;
			rhoR=0.;pR=0.;qVR=0.;
			for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
				if ((*fit).first>=0) { // if contribution is coming from a real cell
					rhoR+=(*fit).second*grid.cell[(*fit).first].rho;
					qVR+=(*fit).second*grid.cell[(*fit).first].v;
					pR+=(*fit).second*grid.cell[(*fit).first].p;
				} else { // if contribution is coming from a ghost cell
					rhoR+=(*fit).second*grid.ghost[-1*((*fit).first+1)].rho;
					qVR+=(*fit).second*grid.ghost[-1*((*fit).first+1)].v;
					pR+=(*fit).second*grid.ghost[-1*((*fit).first+1)].p;
				}
			}
			qVR*=rhoR;

			qNR[0]=rhoR;
			qNR[1]=qVR.dot(grid.face[f].normal);
			qNR[2]=qVR.dot(faceTangent1);
			qNR[3]=qVR.dot(faceTangent2);
			qNR[4]=0.5*(qVR.dot(qVR))/rhoR +pR/ (Gamma - 1.);
			
			//for (int i=0;i<5;++i) qNR[i]=qNL[i]; // outlet condition
			if (bc.region[grid.face[f].bc].type=="outlet" &&
				bc.region[grid.face[f].bc].kind=="fixedPressure") {
 				// find Mach number
				Mach=(qNL[1]/rhoL)/sqrt(Gamma*pL/rhoL);
				if (Mach<1.) pR=bc.region[grid.face[f].bc].p;
				rhoR=rhoL; qVR=qVL;
				qNR[4]=0.5*(qVR.dot(qVR))/rhoR +pR/ (Gamma - 1.);
			}
			if (bc.region[grid.face[f].bc].type=="slip") qNR[1]=-1.*qNL[1];
			if (bc.region[grid.face[f].bc].type=="noslip") {qNR[1]=-1.*qNL[1]; qNR[2]=0.; qNR[3]=0.;}
			if (bc.region[grid.face[f].bc].type=="inlet") {
				rhoR=bc.region[grid.face[f].bc].rho;
				//pR=bc.region[grid.face[f].bc].p;
				pR=pL;
				qVR=rhoR*bc.region[grid.face[f].bc].v;
				qNR[0]=rhoR;
				qNR[1]=qVR.dot(grid.face[f].normal);
				qNR[2]=qVR.dot(faceTangent1);
				qNR[3]=qVR.dot(faceTangent2);
				qNR[4]=0.5*(qVR.dot(qVR))/rhoR +pR/ (Gamma - 1.);
			}

		} else { // partition boundary
			int g=-1*grid.face[f].bc-3;
			rhoR=grid.ghost[g].rho;
			pR=grid.ghost[g].p;
			qVR=rhoR*grid.ghost[g].v;
			qNR[0]=rhoR;
			qNR[1]=qVR.dot(grid.face[f].normal);
			qNR[2]=qVR.dot(faceTangent1);
			qNR[3]=qVR.dot(faceTangent2);
			qNR[4]=0.5*(qVR.dot(qVR))/rhoR +pR/ (Gamma - 1.);
			
		}
		

		roe_flux(qNL, qNR, fluxNormal);

		flux[0] = fluxNormal[0]*grid.face[f].area;
		flux[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
		flux[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
		flux[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
		flux[4] = fluxNormal[4]*grid.face[f].area;

		for (int i=0;i<5;++i) {

			for (int k=0;k<5;++k) temp[k]=qNL[k];
			
			if (i>0 && i<4) {
				deltaV=0.;
				deltaV.comp[i-1]=deltaL[i];
				qNL[1]=(qVL+deltaV).dot(grid.face[f].normal);
				qNL[2]=(qVL+deltaV).dot(faceTangent1);
				qNL[3]=(qVL+deltaV).dot(faceTangent2);
			} else {
				qNL[i]+=deltaL[i];
			}

			roe_flux(qNL, qNR, fluxNormal);

			for (int k=0;k<5;++k) qNL[k]=temp[k];

			fluxPlus[0] = fluxNormal[0]*grid.face[f].area;
			fluxPlus[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
			fluxPlus[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
			fluxPlus[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
			fluxPlus[4] = fluxNormal[4]*grid.face[f].area;

			for (int j=0;j<5;++j) {
				row=grid.cell[parent].globalId*5+j;
				col=grid.cell[parent].globalId*5+i;
				value=(fluxPlus[j]-flux[j])/deltaL[i];
				MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				if (grid.face[f].bc==-1) {
					row=grid.cell[neighbor].globalId*5+j;
					value*=-1.;
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				}
			} // for j


		} // for i

		if (grid.face[f].bc==-1) {

			for (int i=0;i<5;++i) {

				for (int k=0;k<5;++k) temp[k]=qNR[k];
			
				if (i>0 && i<4) {
					deltaV=0.;
					deltaV.comp[i-1]=deltaR[i];
					qNR[1]=(qVR+deltaV).dot(grid.face[f].normal);
					qNR[2]=(qVR+deltaV).dot(faceTangent1);
					qNR[3]=(qVR+deltaV).dot(faceTangent2);
				} else {
					qNR[i]+=deltaR[i];
				}

				roe_flux(qNL, qNR, fluxNormal);

				for (int k=0;k<5;++k) qNR[k]=temp[k];

				fluxPlus[0] = fluxNormal[0]*grid.face[f].area;
				fluxPlus[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
				fluxPlus[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
				fluxPlus[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
				fluxPlus[4] = fluxNormal[4]*grid.face[f].area;

				for (int j=0;j<5;++j) {
					row=grid.cell[neighbor].globalId*5+j;
					col=grid.cell[neighbor].globalId*5+i;
					value=(flux[j]-fluxPlus[j])/deltaR[i];
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
					row=grid.cell[parent].globalId*5+j;
					value*=-1.;
					MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
				} // for j
			} // for i
		} // if internal faces		
		
	} // for faces

	MatAssemblyBegin(impOP,MAT_FLUSH_ASSEMBLY);
	MatAssemblyEnd(impOP,MAT_FLUSH_ASSEMBLY);

	return;
} // end function


