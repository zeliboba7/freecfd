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

void NavierStokes::apply_bcs(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face) {

	bool is_slip=false;
	if (bc[gid][face.bc].kind==SLIP) is_slip=true;

	if (bc[gid][face.bc].type==INLET) {
		if (bc[gid][face.bc].kind==VELOCITY) velocity_inlet(left,right,face);
		else if (bc[gid][face.bc].kind==MDOT) mdot_inlet(left,right,face);
		else if (bc[gid][face.bc].kind==STAGNATION) stagnation_inlet(left,right,face);
	} else if (bc[gid][face.bc].type==OUTLET) {
		outlet(left,right,face);
	} else if (bc[gid][face.bc].type==WALL) {
		wall(left,right,face,is_slip);
	} else if (bc[gid][face.bc].type==SYMMETRY) {
		symmetry(left,right,face);
	}

	return;
} // end apply_bcs

void NavierStokes::velocity_inlet(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face) {
	
	right.V=V.bc(face.bc,face.index);
	// The right cell doesn't exist, exrapolate linearly to an imaginary cell center
	right.V_center=2.*right.V-left.V_center; 
	
	// Note that face normal is pointing out (towards the right state)
	// Extrapolate outgoing Riemann invariant (isentropic) (uN+2a/(gamma-1)
	double uNL=left.V.dot(face.normal);
	double uNR=right.V.dot(face.normal);
	right.a=left.a+0.5*(uNL-uNR)*(material.gamma-1.);
	double MachL=-uNL/left.a;
	if (MachL>=1.) { // supersonic inlet, can't extrapolate anything
		if (bc[gid][face.bc].specified!=BC_STATE) {
			cerr << "[E] Inlet at BC_" << face.bc+1 << " became supersonic, need to specify the themodynamic state" << endl;
			exit(1);
		}
	}
	
	if (bc[gid][face.bc].specified==BC_STATE) {
		right.p=p.bc(face.bc,face.index);
		right.T=T.bc(face.bc,face.index);
		right.T_center=2.*right.T-left.T_center;
		right.rho=rho.bc(face.bc,face.index);
	} else if (bc[gid][face.bc].specified==BC_P) {
		right.p=p.bc(face.bc,face.index);
		right.rho=material.gamma*(right.p+material.Pref)/(right.a*right.a);
		right.T=material.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	} else if (bc[gid][face.bc].specified==BC_T) {
		right.T=T.bc(face.bc,face.index);
		right.T_center=2.*right.T-left.T_center;
		// Extrapolate entropy to get density
		right.rho=left.rho*pow((left.T+material.Tref)/(right.T+material.Tref),1./(1.-material.gamma));
		right.p=material.p(right.rho,right.T);
	} else if (bc[gid][face.bc].specified==BC_RHO) {
		right.rho=rho.bc(face.bc,face.index);
		right.p=right.a*right.a*right.rho/material.gamma-material.Pref;
		right.T=material.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	} else {
		// If nothing is specified extrapolate pressure, get the density
		right.p=left.p;
		right.rho=material.gamma*(right.p+material.Pref)/(right.a*right.a);
		right.T=material.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	}

	return;
} // end velocity_inlet

void NavierStokes::mdot_inlet(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face) {
	
	// Note that face normal is pointing out (towards the right state)
	// Extrapolate outgoing Riemann invariant (isentropic) (uN+2a/(gamma-1)
	double uNL=left.V.dot(face.normal);
	double mdotNR=-1.*(mdot.bc(face.bc,face.index)/bc[gid][face.bc].area)*bc[gid][face.bc].areaVec.norm().dot(face.normal);

	double MachL=-uNL/left.a;
	if (MachL>=1.) { // supersonic inlet, can't extrapolate anything
		if (!bc[gid][face.bc].specified!=BC_STATE) {
			cerr << "[E] Inlet at BC_" << face.bc+1 << " became supersonic, need to specify the themodynamic state" << endl;
			exit(1);
		}
	}
	
	if (bc[gid][face.bc].specified==BC_STATE) {
		right.p=p.bc(face.bc,face.index);
		right.T=T.bc(face.bc,face.index);
		right.T_center=2.*right.T-left.T_center;
		right.rho=rho.bc(face.bc,face.index);
	} else if (bc[gid][face.bc].specified==BC_P) {
		cerr << "[E] You need to specify either T or p with mdot for inlet BC_" << face.bc+1 << endl;
		exit(1);
	} else if (bc[gid][face.bc].specified==BC_T) {
		right.T=T.bc(face.bc,face.index);
		right.T_center=2.*right.T-left.T_center;
		right.a=sqrt(material.gamma*material.R*(right.T+material.Tref));
		if (fabs(mdotNR/uNL-left.rho)/left.rho<0.2) {
			right.rho=mdotNR/(uNL+2.*(left.a-right.a)/(material.gamma-1.));	
		} else {
			right.rho=left.rho;
		}
		right.p=right.a*right.a*right.rho/material.gamma-material.Pref;
	} else if (bc[gid][face.bc].specified==BC_RHO) {
		right.rho=rho.bc(face.bc,face.index);
		right.a=0.5*(material.gamma-1.)*(uNL+2.*left.a/(material.gamma-1.)-mdotNR/right.rho);
		right.p=right.a*right.a*right.rho/material.gamma-material.Pref;
		right.T=material.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
	} else {
		// extrapolate temperature
		right.T=left.T;
		right.T_center=2.*right.T-left.T_center;
		right.a=sqrt(material.gamma*material.R*(right.T+material.Tref));
		if (fabs(mdotNR/uNL-left.rho)/left.rho<0.2) {
			right.rho=mdotNR/(uNL+2.*(left.a-right.a)/(material.gamma-1.));	
		} else {
			right.rho=left.rho;
		}
		right.p=right.a*right.a*right.rho/material.gamma-material.Pref;
	}
	
	right.V=mdotNR/right.rho*face.normal;
	right.V_center=2.*right.V-left.V_center; 

	return;
} // end mdot_inlet

void NavierStokes::stagnation_inlet(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face) {
	// Extrapolate velocity from inside
	right.V=left.V;
	right.V_center=2.*right.V-left.V_center;
	// Impose pressure and temperature
	right.T=T_total.bc(face.bc)-0.5*(material.gamma-1.)/(material.gamma*material.R)*right.V.dot(right.V);

	right.p=(p_total.bc(face.bc)+material.Pref)/(1.+0.5*right.V.dot(right.V)/(material.R*(right.T+material.Tref)))-material.Pref;
	right.T_center=right.T;
	right.rho=material.rho(right.p,right.T);
	//right.T_center=2.*right.T-left.T_center;
	//cout << right.T << "\t" << left.T_center << endl;
	return;
} // end stagnation_inlet

void NavierStokes::outlet(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face) {
	
	double uNL=left.V.dot(face.normal);
	double uNR;
	double MachL=uNL/left.a;
	bool supersonic=false;
	if (MachL>1.) supersonic=true;

	if (bc[gid][face.bc].specified==BC_STATE) {
		cerr << "[E] You can't specify complete thermodynamic state at outlet BC_" << face.bc+1 << endl;
		exit(1);
	} else if (bc[gid][face.bc].specified==BC_P && !supersonic) {
		right.p=p.bc(face.bc,face.index);
		if (MachL<0.) { // reverse flow
			if (bc[gid][face.bc].kind==DAMP_REVERSE) right.p-=0.5*left.rho*uNL*uNL;
		}
		// Extrapolate entropy
		right.rho=left.rho*pow((right.p+material.Pref),1./material.gamma)/pow((left.p+material.Pref),1./material.gamma);
		right.a=sqrt(material.gamma*(right.p+material.Pref)/right.rho);
		// Extrapolate outgoing characteristic
		uNR=uNL+2.*(left.a-right.a)/(material.gamma-1.);
		right.V=left.V-left.V.dot(face.normal)*face.normal+uNR*face.normal;
		//right.V_center=2.*right.V-left.V_center;
		right.V_center=right.V;
		right.T=material.T(right.p,right.rho);
		//right.T_center=2.*right.T-left.T_center;
		right.T_center=right.T;
	} else if (bc[gid][face.bc].specified==BC_T && !supersonic) {
		right.T=T.bc(face.bc,face.index);
		right.T_center=right.T;
		material.gamma=material.gamma; // TODO check left-right gamma's and correct below
		right.a=sqrt(material.gamma*material.R*(right.T+material.Tref));
		// Extrapolate entropy
		right.rho=pow(material.gamma*(left.p+material.Pref)/(right.a*right.a*pow(left.rho,material.gamma)),1.-material.gamma);
		right.p=material.p(right.rho,right.T);
		// Extrapolate outgoing characteristic
		uNR=uNL+2.*(left.a-right.a)/(material.gamma-1.);
		right.V=left.V-left.V.dot(face.normal)*face.normal+uNR*face.normal;
	} else if (bc[gid][face.bc].specified==BC_RHO && !supersonic) {
		right.rho=rho.bc(face.bc,face.index);
		// Extrapolate entropy
		material.gamma=material.gamma; // TODO check left-right gamma's and correct below
		right.p=(left.p+material.Pref)*pow(right.rho/left.rho,material.gamma)-material.Pref;
		right.T=material.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
		right.a=sqrt(material.gamma*(right.p+material.Pref)/right.rho);
		// Extrapolate outgoing characteristic
		uNR=uNL+2.*(left.a-right.a)/(material.gamma-1.);
		right.V=left.V-left.V.dot(face.normal)*face.normal+uNR*face.normal;
	} else {

		// If nothing is specified, everything is extrapolated
		right.p=left.p;
		if (MachL<0.) { // reverse flow
			if (bc[gid][face.bc].kind==DAMP_REVERSE) right.p-=0.5*left.rho*uNL*uNL;
		}
		// Extrapolate entropy
		right.rho=left.rho*pow((right.p+material.Pref),1./material.gamma)/pow((left.p+material.Pref),1./material.gamma);
		right.a=sqrt(material.gamma*(right.p+material.Pref)/right.rho);
		// Extrapolate outgoing characteristic
		uNR=uNL+2.*(left.a-right.a)/(material.gamma-1.);
		right.V=left.V-left.V.dot(face.normal)*face.normal+uNR*face.normal;
		right.V_center=2.*right.V-left.V_center;
		right.T=material.T(right.p,right.rho);
		right.T_center=2.*right.T-left.T_center;
/*		
		right.p=left.p;
		right.T=left.T;
		right.V=left.V;
		right.rho=left.rho;
		right.T_center=right.T;
		right.V_center=right.V;
		//right.V_center=2.*right.V-left.V_center;
		//right.T_center=2.*right.T-left.T_center;
*/
    }
	
	if (MachL<0.) { // reverse flow
		if (bc[gid][face.bc].kind==NO_REVERSE) {
			right.V=-1.*left.V;
			right.V_center=-1.*left.V_center; 
		}
	}
	
	return;
} // end outlet

void NavierStokes::wall(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face,bool is_slip) {
	
	if (bc[gid][face.bc].thermalType==FIXED_T) {
		right.T=T.bc(face.bc,face.index);
		right.T_center=2.*right.T-left.T_center;
		right.p=left.p; // pressure is extrapolated
		right.rho=material.rho(right.p,right.T);
		if (right.rho<=0) {
			right.T=right.T_center;
			right.rho=material.rho(right.p,right.T);
		}
	} else {
		// If nothing is specified, everything is extrapolated
		right.p=left.p;
		right.rho=left.rho;
		right.T=left.T;
		// if temperature is not specified, it is assumed adiabatic
		right.T_center=left.T_center;
	}
	if (is_slip) {
		right.V=left.V-2.*left.V.dot(face.normal)*face.normal;
		right.V_center=left.V_center-2.*left.V_center.dot(face.normal)*face.normal;
	} else {
		right.V=-1.*left.V;
		right.V_center=-1.*left.V_center;
	}

	return;
} // end wall

void NavierStokes::symmetry(NS_Cell_State &left,NS_Cell_State &right,NS_Face_State &face) {

	right.p=left.p;
	right.rho=left.rho;
	right.T=left.T;
	right.T_center=left.T_center;
	
	right.V=left.V-2.*left.V.dot(face.normal)*face.normal;
	right.V_center=left.V_center-2.*left.V_center.dot(face.normal)*face.normal;
	
	return;
} // end symmetry