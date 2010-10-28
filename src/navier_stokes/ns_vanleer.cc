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

void vanLeer_flux(NS_Cell_State &left,NS_Cell_State &right,double fluxNormal[],double Gamma,double Pref,double &weightL) {

	double ML=left.Vn[0]/left.a;
	double MR=right.Vn[0]/right.a;
	double M=0.5*(ML+MR);
	double mdot;
	if (M<=-1.) {
		// Compute fluxes from right
		mdot=right.rho*right.Vn[0];
		fluxNormal[0]=mdot;
		fluxNormal[1]=mdot*right.Vn[0]+right.p;
		fluxNormal[2]=mdot*right.Vn[1];
		fluxNormal[3]=mdot*right.Vn[2];
		fluxNormal[4]=mdot*right.H;
		weightL=0.;
	} else if (M>=1.) {
		// Compute fluxes from left
		mdot=left.rho*left.Vn[0];
		fluxNormal[0]=mdot;
		fluxNormal[1]=mdot*left.Vn[0]+left.p;
		fluxNormal[2]=mdot*left.Vn[1];
		fluxNormal[3]=mdot*left.Vn[2];
		fluxNormal[4]=mdot*left.H;
		weightL=1.;
	} else {
		double fplus=0.25*left.rho*left.a*(ML+1.)*(ML+1.);
		double fminus=-0.25*right.rho*right.a*(1.-MR)*(1.-MR);
		
		fluxNormal[0]=fplus+fminus;
		double termL=(2.*left.a/Gamma)*(0.5*(Gamma-1)*ML+1.);
		double termR=(2.*right.a/Gamma)*(0.5*(Gamma-1)*MR-1.);
		fluxNormal[1]=fplus*termL+fminus*termR-Pref;
		fluxNormal[2]=fplus*left.Vn[1]+fminus*right.Vn[1];
		fluxNormal[3]=fplus*left.Vn[2]+fminus*right.Vn[2];
		termL=0.5*Gamma*Gamma/(Gamma*Gamma-1.)*termL*termL+0.5*(left.Vn[1]*left.Vn[1]+left.Vn[2]*left.Vn[2]);
		termR=0.5*Gamma*Gamma/(Gamma*Gamma-1.)*termR*termR+0.5*(right.Vn[1]*right.Vn[1]+right.Vn[2]*right.Vn[2]);
		fluxNormal[4]=fplus*termL+fminus*termR;
		weightL=0.5;
	}

	return;
} // end vanLeer_flux

