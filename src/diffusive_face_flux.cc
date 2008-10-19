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
#include "inputs.h"
#include "state_cache.h"

extern Grid grid;
extern BC bc;
extern InputFile input;
extern int Rank;
extern double Gamma;
extern double Pref;

void diffusive_face_flux(Cell_State &left,Cell_State &right,Face_State &face,unsigned int f,double flux[]) {

	Vec3D tau_x,tau_y,tau_z,areaVec;
//	double mu;
	
// 	double Cmu=0.09;
// 	double SigmaOmega=1.3;
// 	double SigmaK=1.;
// 	double C1Omega=1.44;
// 	double C2Omega=1.92;
//	double mu_t=0.;




	// Eddy viscosity at the face
	//mu_t=Cmu*faceRho*faceK*faceK/faceomega;
	//mu_t=face.mu;

	areaVec=face.normal*face.area;
	tau_x.comp[0]=2./3.*face.mu* (2.*face.gradU.comp[0]-face.gradV.comp[1]-face.gradW.comp[2]);
	tau_x.comp[1]=face.mu*(face.gradU.comp[1]+face.gradV.comp[0]);
	tau_x.comp[2]=face.mu*(face.gradU.comp[2]+face.gradW.comp[0]);
	tau_y.comp[0]=tau_x.comp[1];
	tau_y.comp[1]=2./3.*face.mu* (2.*face.gradV.comp[1]-face.gradU.comp[0]-face.gradW.comp[2]);
	tau_y.comp[2]=face.mu*(face.gradV.comp[2]+face.gradW.comp[1]);
	tau_z.comp[0]=tau_x.comp[2];
	tau_z.comp[1]=tau_y.comp[2];
	tau_z.comp[2]=2./3.*face.mu*(2.*face.gradW.comp[2]-face.gradU.comp[0]-face.gradV.comp[1]);

	flux[1]+=tau_x.dot(areaVec);
	flux[2]+=tau_y.dot(areaVec);
	flux[3]+=tau_z.dot(areaVec);
	flux[4]+=tau_x.dot(face.v)*areaVec.comp[0]+tau_y.dot(face.v)*areaVec.comp[1]+tau_z.dot(face.v) *areaVec.comp[2];

	// Diffusive k and omega fluxes
// 	flux[5]+=(face.mu+face.mu_t/SigmaK)*gradKf.dot(areaVec);
// 	flux[6]+=(face.mu+face.mu_t/SigmaOmega)*gradOmegaf.dot(areaVec);

	return;
} // end function

