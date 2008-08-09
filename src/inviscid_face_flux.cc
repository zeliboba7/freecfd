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

extern Grid grid;
extern BC bc;
extern int Rank;
extern double Gamma;
extern double Pref;

void update(double dt);

void inviscid_flux_function(unsigned int f, double qL[], double qR[], double flux[]);

void inviscid_face_flux(unsigned int f,string order,string limiter,double flux[], int perturb_state=-1, int perturb_var=0, double perturb_value=0.) {

	double rhoL,rhoR,pL,pR;
	Vec3D vL,vR,faceTangent1,faceTangent2,deltaV;
	double Mach;
	double delta[5];
	unsigned int parent,neighbor;
	double qNL[5], qNR[5], fluxNormal[5];
	
	parent=grid.face[f].parent; neighbor=grid.face[f].neighbor;
	// Take 2 nodes of the face to find a vector in the face plane
	faceTangent1= (grid.face[f].centroid-grid.face[f].node(0));
	faceTangent1/=fabs(faceTangent1);
	// Cross the tangent vector with the normal vector to get the second tangent
	faceTangent2= (grid.face[f].normal).cross(faceTangent1);

	if (order=="second") {
		for (unsigned int i=0;i<5;++i) delta[i]=(grid.face[f].centroid-grid.cell[parent].centroid).dot(grid.cell[parent].limited_grad[i]);
	} else {
		for (unsigned int i=0;i<5;++i) delta[i]=0.;
	}

	if (perturb_state==0) delta[perturb_var]+=perturb_value; // perturb left state

	deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];
	
	// Set left primitive variables
	rhoL=grid.cell[parent].rho+delta[0];
	pL=grid.cell[parent].p+delta[4];
	vL=grid.cell[parent].v+deltaV;

	// Set right primitive variables
	if (grid.face[f].bc==-1) { // means an internal face
		if (order=="second") {
			for (unsigned int i=0;i<5;++i) delta[i]=(grid.face[f].centroid-grid.cell[neighbor].centroid).dot(grid.cell[neighbor].limited_grad[i]);
		}
		if (perturb_state==1) delta[perturb_var]+=perturb_value; // perturb left state
		deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];
		rhoR=grid.cell[neighbor].rho+delta[0];
		pR=grid.cell[neighbor].p+delta[4];
		vR=grid.cell[neighbor].v+deltaV;
		
	} else if (grid.face[f].bc>=0) { // means a boundary face
		// find Mach number
		Mach=(vL.dot(grid.face[f].normal))/sqrt(Gamma*(pL+Pref)/rhoL);
		// Find face averaged variables
		map<int,double>::iterator fit;
		rhoR=0.;pR=0.;vR=0.;
		for (fit=grid.face[f].average.begin();fit!=grid.face[f].average.end();fit++) {
			if ((*fit).first>=0) { // if contribution is coming from a real cell
				rhoR+=(*fit).second*grid.cell[(*fit).first].rho;
				vR+=(*fit).second*grid.cell[(*fit).first].v;
				pR+=(*fit).second*grid.cell[(*fit).first].p;
			} else { // if contribution is coming from a ghost cell
				rhoR+=(*fit).second*grid.ghost[-1*((*fit).first+1)].rho;
				vR+=(*fit).second*grid.ghost[-1*((*fit).first+1)].v;
				pR+=(*fit).second*grid.ghost[-1*((*fit).first+1)].p;
			}
		}
		if (bc.region[grid.face[f].bc].type=="outlet") {
			if (bc.region[grid.face[f].bc].kind=="fixedPressure") {
				if (Mach<1.) pR=bc.region[grid.face[f].bc].p;
			}
		}
		if (bc.region[grid.face[f].bc].type=="slip") {
			rhoR=rhoL;
			vR=vL-2.*vL.dot(grid.face[f].normal)*grid.face[f].normal;
			pR=pL;
		}
		if (bc.region[grid.face[f].bc].type=="noslip") {
			rhoR=rhoL;
			vR=-1.*vL;
			pR=pL;
		}
		if (bc.region[grid.face[f].bc].type=="inlet") {
			rhoR=bc.region[grid.face[f].bc].rho;
			vR=bc.region[grid.face[f].bc].v;
			if (Mach<=-1.) pR=bc.region[grid.face[f].bc].p;
		}
	} else { // partition boundary
		int g=-1*grid.face[f].bc-3;
		if (order=="second") {
			for (unsigned int i=0;i<5;++i) 
			delta[i]=(grid.face[f].centroid-grid.ghost[g].centroid).dot(grid.ghost[g].limited_grad[i]);
		}
		if (perturb_state==1) delta[perturb_var]+=perturb_value; // perturb right state
		deltaV.comp[0]=delta[1]; deltaV.comp[1]=delta[2]; deltaV.comp[2]=delta[3];
		rhoR=grid.ghost[g].rho+delta[0];
		pR=grid.ghost[g].p+delta[4];
		vR=grid.ghost[g].v+deltaV;
	}

	qNL[0]=rhoL;
	qNL[1]=rhoL*vL.dot(grid.face[f].normal);
	qNL[2]=rhoL*vL.dot(faceTangent1);
	qNL[3]=rhoL*vL.dot(faceTangent2);
	qNL[4]=0.5*rhoL*vL.dot(vL)+pL/(Gamma - 1.);

	qNR[0]=rhoR;
	qNR[1]=rhoR*vR.dot(grid.face[f].normal);
	qNR[2]=rhoR*vR.dot(faceTangent1);
	qNR[3]=rhoR*vR.dot(faceTangent2);
	qNR[4]=0.5*rhoR*vR.dot(vR)+pR/(Gamma - 1.);
		
	inviscid_flux_function(f,qNL,qNR,fluxNormal);

	flux[0] = fluxNormal[0]*grid.face[f].area;
	flux[1] = (fluxNormal[1]*grid.face[f].normal.comp[0]+fluxNormal[2]*faceTangent1.comp[0]+fluxNormal[3]*faceTangent2.comp[0]) * grid.face[f].area;
	flux[2] = (fluxNormal[1]*grid.face[f].normal.comp[1]+fluxNormal[2]*faceTangent1.comp[1]+fluxNormal[3]*faceTangent2.comp[1]) * grid.face[f].area;
	flux[3] = (fluxNormal[1]*grid.face[f].normal.comp[2]+fluxNormal[2]*faceTangent1.comp[2]+fluxNormal[3]*faceTangent2.comp[2]) * grid.face[f].area;
	flux[4] = fluxNormal[4]*grid.face[f].area;

	return;
} // end function

