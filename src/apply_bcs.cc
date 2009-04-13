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

void velocity_inlet(Cell_State &left,Cell_State &right,Face_State &face);
void mdot_inlet(Cell_State &left,Cell_State &right,Face_State &face);
void outlet(Cell_State &left,Cell_State &right,Face_State &face);
void noslip(Cell_State &left,Cell_State &right,Face_State &face);
void slip(Cell_State &left,Cell_State &right,Face_State &face);
void symmetry(Cell_State &left,Cell_State &right,Face_State &face);

void apply_bcs(Cell_State &left,Cell_State &right,Face_State &face) {

	if (FLAMELET) {
		Gamma=0.5*(left.gamma+right.gamma);
		eos.R=0.5*(left.R+right.R);
	}
	
	if (bc.region[face.bc].type==INLET) {
		if (bc.region[face.bc].kind==VELOCITY) velocity_inlet(left,right,face);
		else if (bc.region[face.bc].kind==MDOT) mdot_inlet(left,right,face);
	} else if (bc.region[face.bc].type==OUTLET) {
		outlet(left,right,face);
	} else if (bc.region[face.bc].type==NOSLIP) {
		noslip(left,right,face);
	} else if (bc.region[face.bc].type==SLIP) {
		slip(left,right,face);
	} else if (bc.region[face.bc].type==SYMMETRY) {
		symmetry(left,right,face);
	}

	return;
} // end apply_bcs

void velocity_inlet(Cell_State &left,Cell_State &right,Face_State &face) {
	
	right.v=bc.region[face.bc].v;
	right.v_center=2.*right.v-left.v_center; 
	
	// Note that face normal is pointing out (towards the right state)
	// Extrapolate outgoing Riemann invariant (isentropic) (uN+2a/(gamma-1)
	double uNL=left.v.dot(face.normal);
	double uNR=right.v.dot(face.normal);
	right.a=left.a+0.5*(uNL-uNR)*(Gamma-1.);
	double MachL=-uNL/left.a;
	if (MachL>=1.) { // supersonic inlet, can't extrapolate anything
		if (bc.region[face.bc].specified!=BC_STATE) {
			cerr << "[E] Inlet at BC_" << face.bc+1 << " became supersonic, need to specify the themodynamic state" << endl;
			exit(1);
		}
	}
	
	if (bc.region[face.bc].specified==BC_STATE) {
		right.p=bc.region[face.bc].p;
		right.T=bc.region[face.bc].T;
		right.T_center=2.*right.T-left.T_center;
		right.rho=bc.region[face.bc].rho;
	} else if (bc.region[face.bc].specified==BC_P) {
		right.p=bc.region[face.bc].p;
		right.rho=Gamma*(right.p+Pref)/(right.a*right.a);
		right.T=eos.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	} else if (bc.region[face.bc].specified==BC_T) {
		right.T=bc.region[face.bc].T;
		right.T_center=2.*right.T-left.T_center;
		// Extrapolate entropy to get density
		right.rho=left.rho*pow((left.T+Tref)/(right.T+Tref),1./(1.-Gamma));
		right.p=eos.p(right.rho,right.T);
	} else if (bc.region[face.bc].specified==BC_RHO) {
		right.rho=bc.region[face.bc].rho;
		right.p=right.a*right.a*right.rho/Gamma-Pref;
		right.T=eos.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	} else {
		// If nothing is specified extrapolate pressure, get the density
		right.p=left.p;
		right.rho=Gamma*(right.p+Pref)/(right.a*right.a);
		right.T=eos.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	}

	return;
} // end velocity_inlet

void mdot_inlet(Cell_State &left,Cell_State &right,Face_State &face) {
	
	// Note that face normal is pointing out (towards the right state)
	// Extrapolate outgoing Riemann invariant (isentropic) (uN+2a/(gamma-1)
	double uNL=left.v.dot(face.normal);
	double mdotNR=-1.*(bc.region[face.bc].mdot/bc.region[face.bc].area)*bc.region[face.bc].areaVec.norm().dot(face.normal);

	double MachL=-uNL/left.a;
	if (MachL>=1.) { // supersonic inlet, can't extrapolate anything
		if (!bc.region[face.bc].specified!=BC_STATE) {
			cerr << "[E] Inlet at BC_" << face.bc+1 << " became supersonic, need to specify the themodynamic state" << endl;
			exit(1);
		}
	}
	
	if (bc.region[face.bc].specified==BC_STATE) {
		right.p=bc.region[face.bc].p;
		right.T=bc.region[face.bc].T;
		right.T_center=2.*right.T-left.T_center;
		right.rho=bc.region[face.bc].rho;
	} else if (bc.region[face.bc].specified==BC_P) {
		cerr << "[E] You need to specify either T or p with mdot for inlet BC_" << face.bc+1 << endl;
		exit(1);
	} else if (bc.region[face.bc].specified==BC_T) {
		right.T=bc.region[face.bc].T;
		right.T_center=2.*right.T-left.T_center;
		right.a=sqrt(Gamma*eos.R*(right.T+Tref));
		if (fabs(mdotNR/uNL-left.rho)/left.rho<0.2) {
			right.rho=mdotNR/(uNL+2.*(left.a-right.a)/(Gamma-1.));	
		} else {
			right.rho=left.rho;
		}
		right.p=right.a*right.a*right.rho/Gamma-Pref;
	} else if (bc.region[face.bc].specified==BC_RHO) {
		right.rho=bc.region[face.bc].rho;
		right.a=0.5*(Gamma-1.)*(uNL+2.*left.a/(Gamma-1.)-mdotNR/right.rho);
		right.p=right.a*right.a*right.rho/Gamma-Pref;
		right.T=eos.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	} else {
		// extrapolate temperature
		right.T=left.T;
		right.T_center=2.*right.T-left.T_center;
		right.a=sqrt(Gamma*eos.R*(right.T+Tref));
		if (fabs(mdotNR/uNL-left.rho)/left.rho<0.2) {
			right.rho=mdotNR/(uNL+2.*(left.a-right.a)/(Gamma-1.));	
		} else {
			right.rho=left.rho;
		}
		right.p=right.a*right.a*right.rho/Gamma-Pref;
	}
	
	right.v=mdotNR/right.rho*face.normal;
	right.v_center=2.*right.v-left.v_center; 

	return;
} // end mdot_inlet

void outlet(Cell_State &left,Cell_State &right,Face_State &face) {
	
	double uNL=left.v.dot(face.normal);
	double uNR;
	double MachL;
	
	if (bc.region[face.bc].specified==BC_STATE) {
		cerr << "[E] You can't specify complete thermodynamic state at outlet BC_" << face.bc+1 << endl;
		exit(1);
	} else if (bc.region[face.bc].specified==BC_P) {
		right.p=bc.region[face.bc].p;
		MachL=uNL/left.a;
		if (MachL<0.) { // reverse flow
			if (bc.region[face.bc].kind==DAMP_REVERSE) right.p-=0.5*left.rho*uNL*uNL;
		}
		// Extrapolate entropy
		right.rho=left.rho*pow((right.p+Pref),1./right.gamma)/pow((left.p+Pref),1./left.gamma);
		right.a=sqrt(right.gamma*(right.p+Pref)/right.rho);
		// Extrapolate outgoing characteristic
		uNR=uNL+2.*(left.a-right.a)/(Gamma-1.);
		right.v=left.v-left.v.dot(face.normal)*face.normal+uNR*face.normal;
		right.T=eos.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	} else if (bc.region[face.bc].specified==BC_T) {
		right.T=bc.region[face.bc].T;
		right.T_center=right.T;
		Gamma=left.gamma; // TODO check left-right gamma's and correct below
		right.a=sqrt(Gamma*eos.R*(right.T+Tref));
		// Extrapolate entropy
		right.rho=pow(Gamma*(left.p+Pref)/(right.a*right.a*pow(left.rho,Gamma)),1.-Gamma);
		right.p=eos.p(right.rho,right.T);
		// Extrapolate outgoing characteristic
		uNR=uNL+2.*(left.a-right.a)/(Gamma-1.);
		right.v=left.v-left.v.dot(face.normal)*face.normal+uNR*face.normal;
	} else if (bc.region[face.bc].specified==BC_RHO) {
		right.rho=bc.region[face.bc].rho;
		// Extrapolate entropy
		Gamma=left.gamma; // TODO check left-right gamma's and correct below
		right.p=(left.p+Pref)*pow(right.rho/left.rho,Gamma)-Pref;
		right.T=eos.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
		right.a=sqrt(Gamma*(right.p+Pref)/right.rho);
		// Extrapolate outgoing characteristic
		uNR=uNL+2.*(left.a-right.a)/(Gamma-1.);
		right.v=left.v-left.v.dot(face.normal)*face.normal+uNR*face.normal;
	} else {
		// If nothing is specified, everything is extrapolated
		right.p=left.p;
		right.T=left.T; 
		right.T_center=2.*right.T-left.T_center;
		right.rho=left.rho;
		right.v=left.v;
		right.v_center=2.*right.v-left.v_center; 
	}
	
	if (MachL<0.) { // reverse flow
		if (bc.region[face.bc].kind==NO_REVERSE) {
			right.v=-1.*left.v;
			right.v_center=-1.*left.v_center; 
		}
	}
	
	return;
} // end outlet

void noslip(Cell_State &left,Cell_State &right,Face_State &face) {
	
	if (bc.region[face.bc].specified==BC_T) {
		right.T=bc.region[face.bc].T;
		right.T_center=2.*right.T-left.T_center;
		right.p=left.p; // pressure is extrapolated
		right.rho=eos.rho(right.p,right.T); 
	} else {
		// If nothing is specified, everything is extrapolated
		right.p=left.p;
		right.rho=left.rho;
		right.T=left.T;
		// if temperature is not specified, it is assumed adiabatic
		right.T_center=2.*right.T-left.T_center;
	}
	
	right.v=-1.*left.v;
	right.v_center=-1.*left.v_center;

	return;
} // end noslip

void slip(Cell_State &left,Cell_State &right,Face_State &face) {

	if (bc.region[face.bc].specified==BC_T) {
		right.T=bc.region[face.bc].T;
		right.T_center=2.*right.T-left.T_center;
		right.p=left.p; // pressure is extrapolated
		right.rho=eos.rho(right.p,right.T); 
	} else {
		// If nothing is specified, everything is extrapolated
		right.p=left.p;
		right.rho=left.rho;
		right.T=left.T;
		// if temperature is not specified, it is assumed adiabatic
		right.T_center=2.*right.T-left.T_center;
	}
	
	right.v=left.v-2.*left.v.dot(face.normal)*face.normal;
	right.v_center=left.v_center-2.*left.v_center.dot(face.normal)*face.normal;
	
	return;
} // end slip

void symmetry(Cell_State &left,Cell_State &right,Face_State &face) {

	right.p=left.p;
	right.rho=left.rho;
	right.T=left.T;
	right.T_center=left.T_center;
	
	right.v=left.v-2.*left.v.dot(face.normal)*face.normal;
	right.v_center=left.v_center-2.*left.v_center.dot(face.normal)*face.normal;
	
	return;
} // end symmetry