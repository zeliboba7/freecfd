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

double beta=0.125;
double alpha;

double Mach_split_2_plus (double Mach);
double Mach_split_2_minus (double Mach);
double Mach_split_4_plus (double Mach);
double Mach_split_4_minus (double Mach);
double p_split_5_plus (double Mach);
double p_split_5_minus (double Mach);

void AUSMplusUP_flux(NS_Cell_State &left,NS_Cell_State &right,double fluxNormal[],double Gamma,double Pref,double Minf,double &weightL) {

	double Kp=0.25;
	double Ku=0.75;
	double sigma=1.;
	double rho,p,a,M,mdot,Mbar2;
	double aL_hat,aR_hat,aL_star,aR_star;
	double ML,MR;
	double fa=0.;
	double Mref;

	aL_star=left.a;
	aR_star=right.a;
	
	aL_hat=aL_star*aL_star/max(aL_star,left.Vn[0]);
	aR_hat=aR_star*aR_star/max(aR_star,-1.*right.Vn[0]);
	a=min(aL_hat,aR_hat);

	rho=0.5*(left.rho+right.rho);
	ML=left.Vn[0]/a;
	MR=right.Vn[0]/a;

	//Mbar2=0.5*(left.Vn[0]*left.Vn[0]+right.Vn[0]*right.Vn[0])/(a*a);
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
	  -Ku*p_split_5_plus(ML)*p_split_5_minus(MR)*(left.rho+right.rho)*fa*a*(right.Vn[0]-left.Vn[0]);
	
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

