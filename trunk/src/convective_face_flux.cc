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
#include <iostream>
#include <cmath>
#include "inputs.h"
#include "bc.h"
#include "state_cache.h"
using namespace std;

extern BC bc;

double beta=0.125;
double alpha;

void roe_flux(Cell_State &left,Cell_State &right,double fluxNormal[],double &weightL);
void AUSMplusUP_flux(Cell_State &left,Cell_State &right,double fluxNormal[],double &weightL);
void flux_from_right(Cell_State &right,double fluxNormal[]);
double Mach_split_2_plus (double Mach);
double Mach_split_2_minus (double Mach);
double Mach_split_4_plus (double Mach);
double Mach_split_4_minus (double Mach);
double p_split_5_plus (double Mach);
double p_split_5_minus (double Mach);

void convective_face_flux(Cell_State &left,Cell_State &right,Face_State &face,double flux[]) {

	double fluxNormal[5];

	if (face.bc>=0 && bc.region[face.bc].type==INLET) {
		flux_from_right(right,fluxNormal);
		grid.face[face.index].weightL=0.;
	} else if (face.bc>=0 && bc.region[face.bc].type==SYMMETRY) {
		for (int i=0;i<5;++i) fluxNormal[i]=0.;
		fluxNormal[1]=left.p;
		grid.face[face.index].weightL=1.;
	} else if (face.bc>=0 && bc.region[face.bc].type==NOSLIP) {
		roe_flux(left,right,fluxNormal,grid.face[face.index].weightL);
		grid.face[face.index].weightL=1.;
	} else if (CONVECTIVE_FLUX_FUNCTION==ROE) {
		roe_flux(left,right,fluxNormal,grid.face[face.index].weightL);
	} else if (CONVECTIVE_FLUX_FUNCTION==AUSM_PLUS_UP) {
		AUSMplusUP_flux(left,right,fluxNormal,grid.face[face.index].weightL);
	}


	

// 	if (CONVECTIVE_FLUX_FUNCTION==ROE || face.bc>=0) {
// 		roe_flux(left,right,fluxNormal,grid.face[face.index].weightL);
// 	} else if (CONVECTIVE_FLUX_FUNCTION==AUSM_PLUS_UP) {
// 		AUSMplusUP_flux(left,right,fluxNormal,grid.face[face.index].weightL);
// 	}
	
	flux[0] = fluxNormal[0]*face.area;
	flux[1] = (fluxNormal[1]*face.normal[0]+fluxNormal[2]*face.tangent1[0]+fluxNormal[3]*face.tangent2[0])*face.area;
	flux[2] = (fluxNormal[1]*face.normal[1]+fluxNormal[2]*face.tangent1[1]+fluxNormal[3]*face.tangent2[1])*face.area;
	flux[3] = (fluxNormal[1]*face.normal[2]+fluxNormal[2]*face.tangent1[2]+fluxNormal[3]*face.tangent2[2])*face.area;
	flux[4] = fluxNormal[4]*face.area;
	grid.face[face.index].mdot=fluxNormal[0];
// 	if (FLAMELET) flux[4]=0.;

	return;
} // end face flux

void roe_flux(Cell_State &left,Cell_State &right,double fluxNormal[],double &weightL) {

	// Local variables
	double rho,u,v,w,H,a;
	double Du,Dv,Dw,Dp,Drho;
	double mdot;

	double lambda1,lambda5,alpha1,alpha5;
	double right11,right21,right31,right41,right51,right15,right25,right35,right45,right55;
	double gmM1=Gamma-1.;

	// The Roe averaged values
	rho=sqrt(right.rho/left.rho);
	u=(left.vN[0]+rho*right.vN[0])/(1.+rho);
	v=(left.vN[1]+rho*right.vN[1])/(1.+rho);
	w=(left.vN[2]+rho*right.vN[2])/(1.+rho);
	H=(left.H+rho*right.H)/(1.+rho);

	double Denom,GmL,GmR,GH,Gu,Gv,Gw;
	//GmL and GmR are already gamma-1
	Denom=left.H-0.5*left.vN.dot(left.vN);
	GmL=(left.a)/Denom;
	GmR=(right.a)/Denom;
	GH=(GmL*left.H+GmR*rho*right.H)/(1.+rho);
	Gu= (sqrt(GmL)*left.vN[0]+sqrt(GmR)*rho*right.vN[0])/(1.+rho);
	Gv= (sqrt(GmL)*left.vN[1]+sqrt(GmR)*rho*right.vN[1])/(1.+rho);
	Gw= (sqrt(GmL)*left.vN[2]+sqrt(GmR)*rho*right.vN[2])/(1.+rho);
	a=sqrt(GH-0.5*(Gu*Gu+Gv*Gv+Gw*Gw));

	rho=rho*left.rho;

	Du=right.vN[0]-left.vN[0];
	Dv=right.vN[1]-left.vN[1];
	Dw=right.vN[2]-left.vN[2];
	Dp=right.p-left.p;
	Drho=right.rho-left.rho;

	if (u>=0.) { // Calculate from the left side
		lambda1=u-a;
		// Entropy fix
		double Dlambda1=0.5* (min(a,max(0.,2.*(right.a-left.a-Du))));
		if (lambda1<= (-0.5*Dlambda1)) {
			; // do nothing
		} else if (lambda1< (0.5*Dlambda1)) {    // Needs fixin'
			lambda1=-0.5* (lambda1-0.5*Dlambda1)*(lambda1-0.5*Dlambda1) /Dlambda1;
		} else { // just use left values
			lambda1=0.;
		}
		right11=1.;
		right21=u-a;
		right31=v;
		right41=w;
		right51=H-u*a;
		alpha1=0.5*(Dp-rho*a*Du)/(a*a);
		mdot=left.rho*left.vN[0];
		fluxNormal[0]=mdot+lambda1*alpha1*right11;
		fluxNormal[1]=mdot*left.vN[0]+left.p+lambda1*alpha1*right21;
		fluxNormal[2]=mdot*left.vN[1]+lambda1*alpha1*right31;
		fluxNormal[3]=mdot*left.vN[2]+lambda1*alpha1*right41;
		fluxNormal[4]=mdot*left.H+lambda1*alpha1*right51;

	} else { // Calculate from the right side
		lambda5=u+a;
		// Entropy fix
		double Dlambda5=0.5* (min(a,max(0.,2.* (left.a-right.a-Du))));
		if (lambda5<= (-0.5*Dlambda5)) {
			lambda5=0.;
		} else if (lambda5< (0.5*Dlambda5)) {    // Needs fixin'
			lambda5=0.5* (lambda5+0.5*Dlambda5)*(lambda5+0.5*Dlambda5)/Dlambda5;
		} else { // just use left values
			; // do nothing
		}
		right15=1.;
		right25=u+a;
		right35=v;
		right45=w;
		right55=H+u*a;
		alpha5=0.5*(Dp+rho*a*Du)/(a*a);
		
		mdot=right.rho*right.vN[0];
		fluxNormal[0]=mdot-lambda5*alpha5*right15;
		fluxNormal[1]=mdot*right.vN[0]+left.p-lambda5*alpha5*right25;
		fluxNormal[2]=mdot*right.vN[1]-lambda5*alpha5*right35;
		fluxNormal[3]=mdot*right.vN[2]-lambda5*alpha5*right45;
		fluxNormal[4]=mdot*right.H-lambda5*alpha5*right55;
	}
// 	if (abs(left.rho-right.rho)<10e-15) weightL=0.5;
// 	else weightL=((fluxNormal[0]*fluxNormal[0]/fluxNormal[1])-right.rho)/(left.rho-right.rho); // Based on density 
	
// 	if (abs(left.rho*left.vN[0]-right.rho*right.vN[0])<1.e-15) {
// 		weightL=0.5;
// 	} else {
// 		//weightL=(fluxNormal[0]-right.rho*right.vN[0])/(left.rho*left.vN[0]-right.rho*right.vN[0]);
// 		weightL=((fluxNormal[0]*fluxNormal[0]/fluxNormal[1])-right.rho)/(left.rho-right.rho);
// 	}
	weightL=0.;
	
	return;
} // end roe_flux

void AUSMplusUP_flux(Cell_State &left,Cell_State &right,double fluxNormal[],double &weightL) {

	double Kp=0.25;
	double Ku=0.75;
	double sigma=1.;
	double rho,u,p,a,M,mdot,Mbar2,Mo2;
	double aL_hat,aR_hat,Ht,aL_star,aR_star;
	double ML,MR;
	double gmM1=Gamma-1.;
	double fa=0.;
	double Mref;

// 	aL_star=sqrt(2.*gmM1/(Gamma+1.)*left.H);
// 	aR_star=sqrt(2.*gmM1/(Gamma+1.)*right.H);
	aL_star=left.a;
	aR_star=right.a;
	
	aL_hat=aL_star*aL_star/max(aL_star,left.vN[0]);
	aR_hat=aR_star*aR_star/max(aR_star,-1.*right.vN[0]);
	a=min(aL_hat,aR_hat);

	rho=0.5*(left.rho+right.rho);
	ML=left.vN[0]/a;
	MR=right.vN[0]/a;

	//Mbar2=0.5*(left.vN[0]*left.vN[0]+right.vN[0]*right.vN[0])/(a*a);
	Mbar2=0.5*(ML*ML+MR*MR);

 	Mref=max(Minf,sqrt(Mbar2));
 	Mref=min(Mref,1.);

	//Mref=Minf;
	
	if (Mbar2>=1.) {
		fa=1.;
	} else {
		double Mo=sqrt(min(1.,max(Mbar2,Mref*Mref)));
		fa=Mo*(2.-Mo);
	}

	M=Mach_split_4_plus(ML)+Mach_split_4_minus(MR)
	  -Kp/fa*max(1.-sigma*Mbar2,0.)*(right.p-left.p)/(rho*a*a);

	if (M>0) {
		mdot=a*M*left.rho;
	} else {
		mdot=a*M*right.rho;
	}
	
	alpha=3./16.*(-4.+5.*fa*fa);
			
	p=p_split_5_plus(ML)*(left.p+Pref)+p_split_5_minus(MR)*(right.p+Pref)
	  -Ku*p_split_5_plus(ML)*p_split_5_minus(MR)*(left.rho+right.rho)*fa*a*(right.vN[0]-left.vN[0]);
	
	p-=Pref;
	
	fluxNormal[0]=mdot;
	if (mdot>0.) { // Calculate from the left side
		fluxNormal[1]=mdot*left.vN[0]+p;
		fluxNormal[2]=mdot*left.vN[1];
		fluxNormal[3]=mdot*left.vN[2];
		fluxNormal[4]=mdot*left.H;
		weightL=1.;
	} else { // Calculate from the right side
		fluxNormal[1]=mdot*right.vN[0]+p;
		fluxNormal[2]=mdot*right.vN[1];
		fluxNormal[3]=mdot*right.vN[2];
		fluxNormal[4]=mdot*right.H;
		weightL=0.;
	}

	return;
} // end AUSMplusUP_flux


double Mach_split_2_plus (double M) {
	return 0.25*(M+1.)*(M+1.);
}

double Mach_split_2_minus (double M) {
	return -0.25*(M-1.)*(M-1.);
}

double Mach_split_4_plus (double M) {

	if (fabs(M)>=1.) {
		return 0.5*(M+fabs(M));
	} else {
		return Mach_split_2_plus(M)*(1.-16.*beta*Mach_split_2_minus(M));
	}

}

double Mach_split_4_minus (double M) {

	if (fabs(M)>=1.) {
		return 0.5*(M-fabs(M));
	} else {
		return Mach_split_2_minus(M)*(1.+16.*beta*Mach_split_2_plus(M));
	}

}

double p_split_5_plus (double M) {

	if (fabs(M)>=1.) {
		return 0.5*(M+fabs(M))/M;
	} else {
		return Mach_split_2_plus(M)*((2.-M)-16.*alpha*M*Mach_split_2_minus(M));
	}

}

double p_split_5_minus (double M) {

	if (fabs(M)>=1.) {
		return 0.5*(M-fabs(M))/M;
	} else {
		return Mach_split_2_minus(M)*((-2.-M)+16.*alpha*M*Mach_split_2_plus(M));
	}
	
}

void flux_from_right(Cell_State &right,double fluxNormal[]) {

	double mdot=right.rho*right.vN[0];

	fluxNormal[0]=mdot;
	fluxNormal[1]=mdot*right.vN[0]+right.p;
	fluxNormal[2]=mdot*right.vN[1];
	fluxNormal[3]=mdot*right.vN[2];
	fluxNormal[4]=mdot*right.H;

	return;
}
