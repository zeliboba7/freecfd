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

void roe_flux(NS_Cell_State &left,NS_Cell_State &right,double fluxNormal[],double Gamma,double &weightL) {
	
	// Local variables
	double rho,u,v,w,H,a;
	double Du,Dp,Dlambda;
	double lambda,deltaV;
	double mdot,product;
	
	// The Roe averaged values
	rho=sqrt(right.rho/left.rho);
	u=(left.Vn[0]+rho*right.Vn[0])/(1.+rho);
	v=(left.Vn[1]+rho*right.Vn[1])/(1.+rho);
	w=(left.Vn[2]+rho*right.Vn[2])/(1.+rho);
	H=(left.H+rho*right.H)/(1.+rho);
	a=sqrt((Gamma-1.)*(H-0.5*(u*u+v*v+w*w)));
	rho*=left.rho;
	
	Du=right.Vn[0]-left.Vn[0];
	Dp=right.p-left.p;
	
	if (u>=0.) { // Calculate from the left side
		
		deltaV=0.5*(Dp-rho*a*Du)/(a*a); // finite difference of the state vector
		lambda=u-a;
		
		// Entropy fix
		Dlambda=2.*(min(a,max(0.,2.*(left.a-right.a+Du))));
		if (lambda<=(-0.5*Dlambda)) {
			;
		} else if (lambda<(0.5*Dlambda)) {
			lambda=-0.5*(lambda-0.5*Dlambda)*(lambda-0.5*Dlambda)/Dlambda;
		} else { // just use left values
			lambda=0.;
		}
		
		mdot=left.rho*left.Vn[0];
		product=lambda*deltaV;
		
		fluxNormal[0]=mdot+product;
		fluxNormal[1]=mdot*left.Vn[0]+left.p+product*(u-a);
		fluxNormal[2]=mdot*left.Vn[1]+product*v;
		fluxNormal[3]=mdot*left.Vn[2]+product*w;
		fluxNormal[4]=mdot*left.H+product*(H-u*a);

	} else { // Calculate from the right side
		
		deltaV=0.5*(Dp+rho*a*Du)/(a*a); // finite difference of the state vector
		lambda=u+a;
		
		// Entropy fix
		Dlambda=2.*(min(a,max(0.,2.*(right.a-left.a+Du))));
		if (lambda>=(0.5*Dlambda)) {
			;
		} else if (lambda>(-0.5*Dlambda)) {
			lambda=0.5*(lambda+0.5*Dlambda)*(lambda+0.5*Dlambda)/Dlambda;
		} else { // just use right values
			lambda=0.;
		}
		
		mdot=right.rho*right.Vn[0];
		product=lambda*deltaV;
		
		fluxNormal[0]=mdot-product;
		fluxNormal[1]=mdot*right.Vn[0]+right.p-product*(u+a);
		fluxNormal[2]=mdot*right.Vn[1]-product*v;
		fluxNormal[3]=mdot*right.Vn[2]-product*w;
		fluxNormal[4]=mdot*right.H-product*(H+u*a);
	}
	
	// There is straightforward way to determine how stuff is devided for the Roe solver
	weightL=0.5;
	
	return;
} // end roe_flux
