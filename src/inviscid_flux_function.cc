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
#include <iostream>
#include <cmath>
#include <string>
#include "inputs.h"
#include "grid.h"
#include "bc.h"
using namespace std;
extern Grid grid;
extern BC bc;
extern double Gamma;
extern double Pref;
extern string fluxFunction;
extern InputFile input;

double beta=0.125;
double alpha;
double Minf;

void roe_flux(double qL[], double qR[], double flux[]);
void AUSMplusUP_flux(double qL[], double qR[], double flux[]);
double Mach_split_2_plus (double Mach);
double Mach_split_2_minus (double Mach);
double Mach_split_4_plus (double Mach);
double Mach_split_4_minus (double Mach);
double p_split_5_plus (double Mach);
double p_split_5_minus (double Mach);

void inviscid_flux_function(unsigned int f,double qL[],double qR[], double flux[]) {

	if (fluxFunction=="Roe") roe_flux(qL, qR, flux);
	if (fluxFunction=="AUSMplusUP") AUSMplusUP_flux(qL, qR, flux);

	return;
} // end face flux

void roe_flux(double qL[], double qR[], double flux[]) {

	// Local variables
	double rho,u,v,w,h,a;
	double uL,uR,vL,vR,wL,wR,pL,pR,aL,aR,hL,hR;
	double Du,Dv,Dw,Dp,Drho;

	double lambda1,lambda5,alpha1,alpha5;
	double right11,right21,right31,right41,right51,right15,right25,right35,right45,right55;
	double gmM1=Gamma-1.;
	// Assign the left and right velocities to local variables for convenience
	uL=qL[1]/qL[0]; uR=qR[1]/qR[0];
	vL=qL[2]/qL[0]; vR=qR[2]/qR[0];
	wL=qL[3]/qL[0]; wR=qR[3]/qR[0];
	pL=gmM1* (qL[4]-0.5* (qL[1]*uL+qL[2]*vL+qL[3]*wL));
	pR=gmM1* (qR[4]-0.5* (qR[1]*uR+qR[2]*vR+qR[3]*wR));
	aL=sqrt(Gamma*(pL+Pref)/qL[0]); aR=sqrt(Gamma*(pR+Pref)/qR[0]);
	hL=aL*aL/gmM1+0.5* (uL*uL+vL*vL+wL*wL); hR=aR*aR/gmM1+0.5* (uR*uR+vR*vR+wR*wR);
	//pL-=Pref; pR-=Pref;

	// The Roe averaged values
	rho=sqrt(qR[0]/qL[0]);
	u= (uL+rho*uR) / (1.+rho);
	v= (vL+rho*vR) / (1.+rho);
	w= (wL+rho*wR) / (1.+rho);
	h= (hL+rho*hR) / (1.+rho);
	a=sqrt((Gamma-1.) * (h-0.5* (u*u+v*v+w*w)));

	rho=rho*qL[0];

	Du=uR-uL; Dv=vR-vL; Dw=wR-wL; Dp=pR-pL; Drho=qR[0]-qL[0];

	if (u>=0.) { // Calculate from the left side
		lambda1=u-a;
		// Entropy fix
		//double Dlambda1=0.5*(Gamma+1.)*fabs(u-uL);
		double Dlambda1=0.5* (min(a,max(0.,2.* (uL-aL-uR+aR))));
		if (lambda1<= (-0.5*Dlambda1)) {
			; // do nothing
		} else if (lambda1< (0.5*Dlambda1)) {    // Needs fixin'
			lambda1=-0.5* (lambda1-0.5*Dlambda1) * (lambda1-0.5*Dlambda1) /Dlambda1;
		} else { // just use left values
			lambda1=0.;
		}
		right11=1.;
		right21=u-a;
		right31=v;
		right41=w;
		right51=h-u*a;
		alpha1=0.5* (Dp-rho*a*Du) / (a*a);
		flux[0]=qL[0]*uL+lambda1*alpha1*right11;
		flux[1]=qL[1]*uL+pL+lambda1*alpha1*right21;
		flux[2]=qL[1]*vL+lambda1*alpha1*right31;
		flux[3]=qL[1]*wL+lambda1*alpha1*right41;
		flux[4]=uL*(qL[4]+pL+Gamma/(Gamma-1.)*Pref) +lambda1*alpha1*right51;

	} else { // Calculate from the right side
		lambda5=u+a;
		// Entropy fix
		//double Dlambda5=0.5*(Gamma+1.)*fabs(u-uR);
		double Dlambda5=0.5* (min(a,max(0.,2.* (uL+aL-uR-aR))));
		if (lambda5<= (-0.5*Dlambda5)) {
			lambda5=0.;
		} else if (lambda5< (0.5*Dlambda5)) {    // Needs fixin'
			lambda5=0.5* (lambda5+0.5*Dlambda5) * (lambda5+0.5*Dlambda5) /Dlambda5;
		} else { // just use left values
			; // do nothing
		}
		right15=1.;
		right25=u+a;
		right35=v;
		right45=w;
		right55=h+u*a;
		alpha5=0.5* (Dp+rho*a*Du) / (a*a);
		flux[0]=qR[0]*uR-lambda5*alpha5*right15;
		flux[1]=qR[1]*uR+pR-lambda5*alpha5*right25;
		flux[2]=qR[1]*vR-lambda5*alpha5*right35;
		flux[3]=qR[1]*wR-lambda5*alpha5*right45;
		flux[4]=uR*(qR[4]+pR+Gamma/(Gamma-1.)*Pref)-lambda5*alpha5*right55;
	}

	flux[5]=0.; flux[6]=0.; // TODO need to look how to do turbulence convective fluxes in Roe solver

	return;
} // end roe_flux

void AUSMplusUP_flux(double qL[], double qR[], double flux[]) {

	double Kp=0.25;
	double Ku=0.75;
	double sigma=1.;
	// Local variables
	double rho,u,p,a,M,mdot,Mbar2,Mo2;
	double aL_hat,aR_hat,Ht,aL_star,aR_star;
	double rhoL,rhoR,uL,uR,vL,vR,wL,wR,pL,pR,aL,aR,aL2,aR2,ML,MR,HL,HR;
	double kL,kR,omegaL,omegaR;
	double gmM1=Gamma-1.;
	double fa=0.;

	Minf=input.section["numericalOptions"].doubles["Minf"];
	// Assign the left and right velocities to local variables for convenience
	rhoL=qL[0]; rhoR=qR[0];
	uL=qL[1]/rhoL; uR=qR[1]/rhoR;
	vL=qL[2]/rhoL; vR=qR[2]/rhoR;
	wL=qL[3]/rhoL; wR=qR[3]/rhoR;
	pL=gmM1*(qL[4]-0.5*(qL[1]*uL+qL[2]*vL+qL[3]*wL));
	pR=gmM1*(qR[4]-0.5*(qR[1]*uR+qR[2]*vR+qR[3]*wR));
	kL=qL[5]/qL[0]; kR=qR[5]/qR[0];
	omegaL=qL[6]/qL[0]; omegaR=qR[6]/qR[0];
	aL2=Gamma*(pL+Pref)/rhoL; aR2=Gamma*(pR+Pref)/rhoR;
	HL=aL2/gmM1+0.5*(uL*uL+vL*vL+wL*wL); HR=aR2/gmM1+0.5*(uR*uR+vR*vR+wR*wR);
	aL=sqrt(aL2); aR=sqrt(aR2);

	aL_star=sqrt(2.*gmM1/(Gamma+1.)*HL);
	aR_star=sqrt(2.*gmM1/(Gamma+1.)*HR);
	aL_hat=aL_star*aL_star/max(aL_star,uL);
	aR_hat=aR_star*aR_star/max(aR_star,-1.*uR);
	a=min(aL_hat,aR_hat);
	//a=0.5*(aL+aR);
	
	rho=0.5*(rhoL+rhoR);
	ML=uL/a;
	MR=uR/a;

	Mbar2=0.5*(uL*uL+uR*uR)/(a*a);

	if (Mbar2>=1.) {
		fa=1.;
	} else {
		double Mo=sqrt(min(1.,max(Mbar2,Minf*Minf)));
		fa=Mo*(2.-Mo);
	}

	M=Mach_split_4_plus(ML)+Mach_split_4_minus(MR)
	  -Kp/fa*max(1.-sigma*Mbar2,0.)*(pR-pL)/(rho*a*a);
	
	if (M>0) {
		mdot=a*M*rhoL;
	} else {
		mdot=a*M*rhoR;
	}

	alpha=3./16.*(-4.+5.*fa*fa);
			
	p=p_split_5_plus(ML)*pL+p_split_5_minus(MR)*pR
	  -Ku*p_split_5_plus(ML)*p_split_5_minus(MR)*(rhoL+rhoR)*fa*a*(uR-uL);
			
	flux[0]=mdot;
	if (mdot>0.) { // Calculate from the left side
		flux[1]=mdot*uL+p;
		flux[2]=mdot*vL;
		flux[3]=mdot*wL;
		flux[4]=mdot*HL;
		flux[5]=mdot*kL;
		flux[6]=mdot*omegaL;
	} else { //if (M<0.) { // Calculate from the right side
		flux[1]=mdot*uR+p;
		flux[2]=mdot*vR;
		flux[3]=mdot*wR;
		flux[4]=mdot*HR;
		flux[5]=mdot*kR;
		flux[6]=mdot*omegaR;
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
