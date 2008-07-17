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
using namespace std;
extern double Gamma;

void roe_flux(const double qL[], const double qR[], double flux[]) {
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
	aL=sqrt(Gamma*pL/qL[0]); aR=sqrt(Gamma*pR/qR[0]);
	hL=aL*aL/gmM1+0.5* (uL*uL+vL*vL+wL*wL); hR=aR*aR/gmM1+0.5* (uR*uR+vR*vR+wR*wR);

	// The Roe averaged values
	rho=sqrt(qR[0]/qL[0]);
	u= (uL+rho*uR) / (1.+rho);
	v= (vL+rho*vR) / (1.+rho);
	w= (wL+rho*wR) / (1.+rho);
	h= (hL+rho*hR) / (1.+rho);
	a=sqrt((Gamma-1.) * (h-0.5* (u*u+v*v+w*w)));

	rho=rho*qL[0];

	Du=uR-uL; Dv=vR-vL; Dw=wR-wL; Dp=pR-pL; Drho=qR[0]-qL[0];

	if (u>=0) { // Calculate from the left side
		lambda1=u-a;
		// Entropy fix
		//double Dlambda1=0.5*(Gamma+1.)*fabs(u-uL);
		double Dlambda1=0.5* (min(a,max(0.,2* (uL-aL-uR+aR))));
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
		flux[4]=uL* (qL[4]+pL) +lambda1*alpha1*right51;
	} else { // Calculate from the right side
		lambda5=u+a;
		// Entropy fix
		//double Dlambda5=0.5*(Gamma+1.)*fabs(u-uR);
		double Dlambda5=0.5* (min(a,max(0.,2* (uL+aL-uR-aR))));
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
		flux[4]=uR* (qR[4]+pR)-lambda5*alpha5*right55;
	}

	return;
}
