/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

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
#include "petsc_functions.h"
#include "bc.h"
#include "rans.h"

extern BC bc;
extern RANS rans;
		
inline void preconditioner_ws95(Cell &c,int cid,double &mu,double P[][5]) {
	
	double p=c.p+Pref;
	double T=c.T+Tref;
	// Reference velocity
	double Ur;
	
	double c_p=Gamma/(Gamma-1.)*p/(c.rho*T);
	
	// Speed of sound 
	double a=sqrt(Gamma*p/c.rho); // For ideal gas
	
	Ur=max(fabs(c.v),1.e-5*a);
	Ur=min(Ur,a);
	if (EQUATIONS==NS) Ur=max(Ur,mu/(c.rho*c.lengthScale));
	
	double deltaPmax=0.;
	// Now loop through neighbor cells to check pressure differences
	for (it=c.neighborCells.begin();it!=c.neighborCells.end();it++) {
		deltaPmax=max(deltaPmax,fabs(c.p-grid.cell[*it].p));
	}
	// TODO loop ghosts too
	
	Ur=max(Ur,1.e-5*sqrt(deltaPmax/c.rho));
	
	double drho_dp; // Derivative of density w.r.t press. @ const temp.
	double drho_dT; // Derivative of density w.r.t temp. @ const. press
	
	// For ideal gas
	drho_dT=-1.*c.rho/T;

	drho_dp=1./(Ur*Ur)-drho_dT/(c.rho*c_p); // This is the only change over non-preconditioned
	
	double H=c_p*T+0.5*c.v.dot(c.v);
	
	// Conservative to primite Jacobian
	P[0][0]=drho_dp; P[0][4]=drho_dT;
			
	P[1][0]=drho_dp*c.v[0]; P[1][1]=c.rho; P[1][4]=drho_dT*c.v[0];
	P[2][0]=drho_dp*c.v[1]; P[2][2]=c.rho; P[2][4]=drho_dT*c.v[1];
	P[3][0]=drho_dp*c.v[2]; P[3][3]=c.rho; P[3][4]=drho_dT*c.v[2];
		
	P[4][0]=drho_dp*H-1.; 
	P[4][1]=c.rho*c.v[0];
	P[4][2]=c.rho*c.v[1];
	P[4][3]=c.rho*c.v[2];
	P[4][4]=drho_dT*H+c.rho*c_p; 
	
	return;
	
}

void pseudo_time_terms() {
	
	// At this point in the first pseudo time step,
	// both the conservative to primitive conversion
	// and flux jacobians are added to impOP
	
	// Add the contribution to rhs
	if (ps_timeStep==1) {
		VecSet(pseudo_right,0.);
		MatZeroEntries(pseudo_time);
	} else {
		MatMultAdd(impOP,deltaU,pseudo_right,pseudo_right); // pseudo_right+=impOP*deltaU
	}
		
	// Add the unsteady term

	double P [5][5]; // preconditioner
	for (int i=0;i<5;++i) for (int j=0;j<5;++j) P[i][j]=0.;

	PetscInt row,col;
	PetscScalar value;
	
	for (int c=0;c<grid.cellCount;++c) {
		// Only one preconditioner is available for now
		double mu=viscosity;
		preconditioner_ws95(grid.cell[c],c,mu,P);

		for (int i=0;i<5;++i) {
			row=(grid.myOffset+c)*5+i;
			for (int j=0;j<5;++j) {
				col=(grid.myOffset+c)*5+j;
				value=P[i][j]*grid.cell[c].volume/grid.cell[c].dt;
				MatSetValues(pseudo_time,1,&row,1,&col,&value,INSERT_VALUES);
			}
		}
		
	}
	
	return;
}

