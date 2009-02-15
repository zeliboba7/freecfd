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
#include "state_cache.h"

void sources(Cell_State &state,double source[],bool forJacobian=false) {
	if (TURBULENCE_MODEL!=NONE) {
		if (!forJacobian) {
			double alpha=5./9.;
			double beta=3./40.;
			double betaStar=0.09;	
			double mu_t=fabs(state.rho*state.k/state.omega);
			double divU=0.;
			double tauGradU=0.;
			
			divU=state.gradU[0]+state.gradV[1]+state.gradW[2];
			
			tauGradU+=state.gradU.dot(state.gradU)+state.gradV.dot(state.gradV)+state.gradW.dot(state.gradW);
			tauGradU+=state.gradU[0]*state.gradU[0]+state.gradV[1]*state.gradV[1]+state.gradW[2]*state.gradW[2];
			tauGradU+=2.*(state.gradV[0]*state.gradU[1]+state.gradW[0]*state.gradU[2]+state.gradW[1]*state.gradV[2]);
			tauGradU*=mu_t;
			tauGradU-=2./3.*divU*(mu_t*divU+state.rho*state.k);
			
			source[5]=(tauGradU-state.rho*betaStar*state.omega*state.k)*state.volume;
			source[6]=(alpha*state.omega/state.k*tauGradU-state.rho*beta*state.omega*state.omega)*state.volume;
		}
	}	

	return;
} // end function

