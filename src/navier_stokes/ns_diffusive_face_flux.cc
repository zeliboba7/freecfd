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
#include "ns.h"
#include "rans.h"
extern vector<RANS> rans;

void NavierStokes::diffusive_face_flux(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face,double flux[]) {

	Vec3D tau_x,tau_y,tau_z,areaVec;

	double qq;
	double turb_visc=0.;
	double turb_cond=0.;
	
	if (turbulent[gid]) {
		turb_visc=rans[gid].mu_t.face(face.index);
		turb_cond=material.Cp(face.T)*turb_visc/material.Pr;
	}
	
	areaVec=face.normal*face.area;
	tau_x[0]=2./3.*(2.*face.gradu[0]-face.gradv[1]-face.gradw[2]);
	tau_x[1]=face.gradu[1]+face.gradv[0];
	tau_x[2]=face.gradu[2]+face.gradw[0];
	tau_y[0]=tau_x[1];
	tau_y[1]=2./3.* (2.*face.gradv[1]-face.gradu[0]-face.gradw[2]);
	tau_y[2]=face.gradv[2]+face.gradw[1];
	tau_z[0]=tau_x[2];
	tau_z[1]=tau_y[2];
	tau_z[2]=2./3.*(2.*face.gradw[2]-face.gradu[0]-face.gradv[1]);

	flux[1]=(face.mu+turb_visc)*tau_x.dot(areaVec);
	flux[2]=(face.mu+turb_visc)*tau_y.dot(areaVec);
	flux[3]=(face.mu+turb_visc)*tau_z.dot(areaVec);
	flux[4]=(face.mu+turb_visc)*(tau_x.dot(face.V)*areaVec[0]+tau_y.dot(face.V)*areaVec[1]+tau_z.dot(face.V)*areaVec[2]);
	
	qq=(face.lambda+turb_cond)*face.gradT.dot(areaVec);
	if (face.bc>=0) {
		if (bc[gid][face.bc].thermalType==FIXED_Q) qq=qdot.bc(face.bc,face.index)*face.area;
		else if (qdot.bcValue[face.bc].size()>1) qdot.bc(face.bc,face.index)=-qq;
	}
	flux[4]+=qq;
	
	return;
}

