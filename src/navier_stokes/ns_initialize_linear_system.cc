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

// Conservative to primitive variable set conversion matrix
void NavierStokes::cons2prim(int c,vector<vector<double> > &P) {

	double pressure=p.cell(c)+material.Pref;
	double temperature=T.cell(c)+material.Tref;
	double drho_dT; // Derivative of density w.r.t temp. @ const. press
	double drho_dp; // Derivative of density w.r.t press. @ const temp.
	
	// For ideal gas
	drho_dT=-1.*rho.cell(c)/temperature;
	drho_dp=rho.cell(c)/pressure;

	double c_p=material.gamma/(material.gamma-1.)*pressure/(rho.cell(c)*temperature);
	
	double H=c_p*temperature+0.5*V.cell(c).dot(V.cell(c));
	
	// Conservative to primite Jacobian
	P[0][0]=drho_dp; P[0][4]=drho_dT;
	
	P[1][0]=drho_dp*V.cell(c)[0]; P[1][1]=rho.cell(c); P[1][4]=drho_dT*V.cell(c)[0];
	P[2][0]=drho_dp*V.cell(c)[1]; P[2][2]=rho.cell(c); P[2][4]=drho_dT*V.cell(c)[1];
	P[3][0]=drho_dp*V.cell(c)[2]; P[3][3]=rho.cell(c); P[3][4]=drho_dT*V.cell(c)[2];
	
	P[4][0]=drho_dp*H-1.;
	P[4][1]=rho.cell(c)*V.cell(c)[0];
	P[4][2]=rho.cell(c)*V.cell(c)[1];
	P[4][3]=rho.cell(c)*V.cell(c)[2];
	P[4][4]=drho_dT*H+rho.cell(c)*c_p;

	return;
}

void NavierStokes::initialize_linear_system() {

	MatZeroEntries(impOP);

	// Preconditioner, in this instance it is simply
	// a converter from conservative set of variables to primitive 
	vector<vector<double> > P;
	P.resize(5);
	for (int i=0;i<5;++i) {
		P[i].resize(5);
		for (int j=0;j<5;++j) P[i][j]=0.;
	}

	PetscInt row,col;
	PetscScalar value;
	
	for (int c=0;c<grid[gid].cellCount;++c) {

		cons2prim(c,P);

		for (int i=0;i<5;++i) {
			row=(grid[gid].myOffset+c)*5+i;
			for (int j=0;j<5;++j) {
				col=(grid[gid].myOffset+c)*5+j;
				value=P[i][j]*grid[gid].cell[c].volume/dt[gid].cell(c);
				MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			}
		}
		
	}
	
	return;
}

