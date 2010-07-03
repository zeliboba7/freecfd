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

template<typename T>
inline T signum(T n)
{
	if (n < 0) return -1;
	if (n > 0) return 1;
	return 0;
}

double Csd1=0.1;
double Csd2=10.;
double fa=0.;

void SD_SLAU_flux(NS_Cell_State &left,NS_Cell_State &right,double fluxNormal[],double Pref,double &weightL) {

	double a=0.5*(left.a+right.a);
	double Mhat=min(1.,sqrt(0.5*(left.V.dot(left.V)+right.V.dot(right.V))/a));
	double chi=(1.-Mhat)*(1.-Mhat);
	double Mplus=left.Vn[0]/a;
	double Mminus=right.Vn[0]/a;
	double Bplus,Bminus;
	double p,ave_p,delta_p,max_delta_p;
	double alpha=3./16.*(-4.+5.*fa*fa);
	alpha=0.;
	if (fabs(Mplus)<1) {
		Bplus=0.25*(2.-Mplus)*(Mplus+1.)*(Mplus+1.)+alpha*Mplus*(Mplus-1.)*(Mplus-1.);
	} else {
		Bplus=0.5*(1.+signum(Mplus));
	}
	
	if (fabs(Mminus)<1) {
		Bminus=0.25*(2.+Mminus)*(Mminus-1.)*(Mminus-1.)-alpha*Mminus*(Mminus-1.)*(Mminus-1.);
	} else {
		Bminus=0.5*(1.+signum(-Mminus));
	}
	
	ave_p=0.5*(left.p+right.p)+Pref;
	delta_p=right.p-left.p;
	p=ave_p-0.5*(Bplus-Bminus)*delta_p+(1.-chi)*(Bplus+Bminus-1.)*ave_p;
	
	// TODO max_delta_p should be found out from surrounding cells
	max_delta_p=fabs(delta_p);
	
	double g=-max(min(Mplus,0.),-1.)*min(max(Mminus,0.),1.);
	double Vn=(left.rho*fabs(left.Vn[0])+right.rho*fabs(right.Vn[0]))/(left.rho+right.rho);
	double Vnplus=(1.-g)*Vn+g*fabs(left.Vn[0]);
	double Vnminus=(1.-g)*Vn+g*fabs(right.Vn[0]);
	double theta=(Csd2*fabs(delta_p)/ave_p+Csd1)/(max_delta_p/ave_p+Csd1);
	theta=min(1.,theta*theta);
	theta=1.; // See the above TODO to enable shock detection
	double mdot=0.5*(left.rho*(left.Vn[0]+Vnplus)+right.rho*(right.Vn[0]-Vnminus)-theta*max(0.,(1.-Vn/a))*delta_p/a);
	
	p-=Pref;
	
	fluxNormal[0]=mdot;
	if (mdot>0.) { // Calculate from the left side
		fluxNormal[1]=mdot*left.Vn[0]+p;
		fluxNormal[2]=mdot*left.Vn[1];
		fluxNormal[3]=mdot*left.Vn[2];
		fluxNormal[4]=mdot*left.H;
		weightL=1.;
	} else { // Calculate from the right side
		fluxNormal[1]=mdot*right.Vn[0]+p;
		fluxNormal[2]=mdot*right.Vn[1];
		fluxNormal[3]=mdot*right.Vn[2];
		fluxNormal[4]=mdot*right.H;
		weightL=0.;
	}
		
	return;
} // end SD_SLAU_flux
