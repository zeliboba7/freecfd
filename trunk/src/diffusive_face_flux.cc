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
#include "commons.h"
#include <cmath>
#include "bc.h"
#include "state_cache.h"

extern BC bc;

void diffusive_face_flux(Cell_State &left,Cell_State &right,Face_State &face,double flux[]) {

	Vec3D tau_x,tau_y,tau_z,areaVec;
	
	/*
 	double alpha=5./9.;
	double SigmaOmega=0.5;
 	double SigmaK=0.5;
	double mu_t=face.rho*face.k/face.omega;
	mu_t=0.;
	*/
	
	areaVec=face.normal*face.area;
	tau_x[0]=2./3.*(2.*face.gradU[0]-face.gradV[1]-face.gradW[2]);
	tau_x[1]=face.gradU[1]+face.gradV[0];
	tau_x[2]=face.gradU[2]+face.gradW[0];
	tau_y[0]=tau_x[1];
	tau_y[1]=2./3.* (2.*face.gradV[1]-face.gradU[0]-face.gradW[2]);
	tau_y[2]=face.gradV[2]+face.gradW[1];
	tau_z[0]=tau_x[2];
	tau_z[1]=tau_y[2];
	tau_z[2]=2./3.*(2.*face.gradW[2]-face.gradU[0]-face.gradV[1]);

	flux[1]=(viscosity)*tau_x.dot(areaVec);
	flux[2]=(viscosity)*tau_y.dot(areaVec);
	flux[3]=(viscosity)*tau_z.dot(areaVec);
	flux[4]=(viscosity)*(tau_x.dot(face.v)*areaVec[0]+tau_y.dot(face.v)*areaVec[1]+tau_z.dot(face.v)*areaVec[2]);
	flux[4]+=conductivity*face.T*face.area;
	flux[5]=0.;
	flux[6]=0.;

	/*
	// Diffusive k and omega fluxes
 	flux[5]=(viscosity+mu_t/SigmaK)*face.gradK.dot(areaVec);
 	flux[6]=(viscosity+mu_t/SigmaOmega)*face.gradOmega.dot(areaVec);

	double tauUgrad=face.rho/(viscosity+mu_t)*(face.v[0]*flux[1]+face.v[1]*flux[2]+face.v[2]*flux[3]);
	flux[5]=tauUgrad;
	flux[6]=alpha*face.omega/face.k*tauUgrad;
	*/

	return;
} // end function

