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
#include "petsc_functions.h"
#include "bc.h"

extern BC bc;

inline void preconditioner_none(Cell &c, double P[][7]) {
	
	// Conservative to primite Jacobian
	P[0][0]=1.;
	P[1][0]=c.v[0]; P[1][1]=c.rho;
	P[2][0]=c.v[1]; P[2][2]=P[1][1];
	P[3][0]=c.v[2]; P[3][3]=P[1][1];
		
	P[4][0]=0.5*c.v.dot(c.v);
	P[4][1]=P[1][1]*P[1][0];
	P[4][2]=P[1][1]*P[2][0];
	P[4][3]=P[1][1]*P[3][0];
	P[4][4]=1./(Gamma-1.);

	P[5][0]=c.k; P[5][5]=P[1][1];
	P[6][0]=c.omega; P[6][6]=P[1][1];
	
	return;
}
		
// void preconditioner_cm91(unsigned int c, double P[][7]);
void preconditioner_ws95(unsigned int c, double P[][7]);

void mat_print(double P[][7]);

void initialize_linear_system() {

	int nSolVar=5; // Basic equations to solve

	if (TURBULENCE_MODEL!=NONE) nSolVar+=2;

	if ((timeStep) % jacobianUpdateFreq == 0 | timeStep==restart+1) MatZeroEntries(impOP);
	
	PetscInt counter=0;

	double d,lengthScale,dtLocal,a;
	double P [7][7]; // preconditioner
	for (int i=0;i<7;++i) for (int j=0;j<7;++j) P[i][j]=0.;

	PetscInt row,col;
	PetscScalar value;
			
	for (unsigned int c=0;c<grid.cellCount;++c) {

		if (PRECONDITIONER==WS95) {
			preconditioner_ws95(c,P);
		} else {
			preconditioner_none(grid.cell[c],P);
		}

		if (TIME_STEP_TYPE==CFL_LOCAL) {
			// Determine time step with CFL condition
			lengthScale;
			dtLocal=1.E20;
			a=sqrt(Gamma*(grid.cell[c].p+Pref)/grid.cell[c].rho);
			lengthScale=grid.cell[c].lengthScale;
			dtLocal=min(dtLocal,CFLlocal*lengthScale/(fabs(grid.cell[c].v[0])+a));
			dtLocal=min(dtLocal,CFLlocal*lengthScale/(fabs(grid.cell[c].v[1])+a));
			dtLocal=min(dtLocal,CFLlocal*lengthScale/(fabs(grid.cell[c].v[2])+a));
			d=grid.cell[c].volume/dtLocal;
		} else {
			d=grid.cell[c].volume/dt;
		}
		
		for (int i=0;i<nSolVar;++i) {
			row=grid.cell[c].globalId*nSolVar+i;
			for (int j=0;j<nSolVar;++j) {
				col=grid.cell[c].globalId*nSolVar+j;
				value=P[i][j]*d;
				MatSetValues(impOP,1,&row,1,&col,&value,ADD_VALUES);
			}
		}
			
	}

	return;
}


void preconditioner_ws95(unsigned int c, double P[][7]) {

	double Umax,Ur,theta,vis;
	double rho,u,v,w,p,H,q2,a2,q,a,k,omega,rhoT,R,cp;
	double temp [5][5];
				
	rho=grid.cell[c].rho;
	u=grid.cell[c].v.comp[0];
	v=grid.cell[c].v.comp[1];
	w=grid.cell[c].v.comp[2];
	p=grid.cell[c].p+Pref;
	k=grid.cell[c].k;
	omega=grid.cell[c].omega;
	q2=grid.cell[c].v.dot(grid.cell[c].v);
	a2=Gamma*p/rho;
	H=0.5*q2+a2/(Gamma - 1.);
	
	for (int i=0;i<7;++i) for (int j=0;j<7;++j) P[i][j]=0.;
	
	P[0][0]=1.;
	P[0][4]=-1.*rho/p;

	P[1][0]=u;
	P[1][1]=rho;
	P[1][4]=-1.*rho*u/p;
	P[2][0]=v;
	P[2][2]=rho;
	P[2][4]=-1.*rho*v/p;
	P[3][0]=w;
	P[3][3]=rho;
	P[3][4]=-1.*rho*w/p;
	
	P[4][0]=0.5*q2;
	P[4][1]=rho*u;
	P[4][2]=rho*v;
	P[4][3]=rho*w;
	P[4][4]=-1.*rho*H/p;

	P[5][0]=k; P[5][5]=rho;
	P[6][0]=omega; P[6][6]=rho;
	
	return;
}

// void preconditioner_cm91(unsigned int c, double P[][7]) {
// 
// 	double Mach,Mach2,a2;
// 	double rho,u,v,w,p,k,omega,e,d,q2;
// 	double beta;
// 	double gm1=Gamma-1.;
// 	double R=287.;
// 	
// 	rho=grid.cell[c].rho;
// 	u=grid.cell[c].v.comp[0];
// 	v=grid.cell[c].v.comp[1];
// 	w=grid.cell[c].v.comp[2];
// 	p=grid.cell[c].p;
// 	k=grid.cell[c].k;
// 	omega=grid.cell[c].omega;
// 	q2=grid.cell[c].v.dot(grid.cell[c].v);
// 	a2=Gamma*(p+Pref)/rho;
// 	Mach2=q2/a2;
// 	Mach=sqrt(Mach2);
// 	
// 	if (Mach<=1.e-5) {
// 		Mach=1.e-5;
// 	} else if (Mach<1.) {
// 		//use local
// 	} else {
// 		Mach=1.;
// 	}
// 
// 	Mach2=Mach*Mach;
// 
// 	beta=max(1.e-5,q2);
// 			
// 	//beta=a2;
// 	
// 	e=0.5*q2+a2/(Gamma*(Gamma - 1.));
// 	
// 	for (int i=0;i<7;++i) for (int j=0;j<7;++j) P[i][j]=0.;
// 
// 	P[0][0]=1./(beta*Mach2);
// 	
// 	P[1][0]=u/(beta*Mach2);
// 	P[1][1]=rho;
// 
// 	P[2][0]=v/(beta*Mach2);
// 	P[2][2]=rho;
// 
// 	P[3][0]=w/(beta*Mach2);
// 	P[3][3]=rho;
// 
// 	P[4][0]=(rho*e+p)/(rho*beta*Mach2)-1.;
// 	P[4][1]=rho*u;
// 	P[4][2]=rho*v;
// 	P[4][3]=rho*w;
// 	P[4][4]=Gamma*rho*R/(Gamma-1.);
// 
// 	P[5][0]=k; P[5][5]=rho;
// 	P[6][0]=omega; P[6][6]=rho;
// 	
// 	return;
// }

void mat_print(double mat[][5]) {


	cout << endl;
	for (int i=0;i<5;++i) {
		for (int j=0;j<5;++j) {
			cout << mat[i][j] << "\t";
		}
		cout << endl;
	}
	
	return;
}
