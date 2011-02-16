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
#include "material.h"

MATERIAL::MATERIAL() {
	;
}

void MATERIAL::set(int gid) {
	if (input.section("grid",gid).get_string("material").is_found) {
		Mw=material_input[gid].get_double("molarmass");
		density=material_input[gid].get_double("density"); // for solids
		gamma=material_input[gid].get_double("gamma");
		
		if (material_input[gid].section("equationofstate").get_string("model")=="idealgas") {
			eos_model=IDEAL_GAS;
			R=UNIV_GAS_CONST/Mw;
		} else {
			if (Rank==0) cerr << "[grid=" << gid+1 << "] invalid equation of state option!" << endl;
			MPI_Abort(MPI_COMM_WORLD,-1);
		}
		
		if (material_input[gid].get_double("viscosity").is_found) {
			visc_model=CONSTANT;
			mu=material_input[gid].get_double("viscosity");
		} else if (material_input[gid].section("viscosity").is_found) {
			if (material_input[gid].section("viscosity").get_string("model")=="sutherlands") {
				visc_model=SUTHERLANDS;
				sut_mu_ref=material_input[gid].section("viscosity").get_double("referenceviscosity");
				sut_T_ref=material_input[gid].section("viscosity").get_double("referencetemperature");
				sut_S=material_input[gid].section("viscosity").get_double("sutherlandtemperature");
			} else {
				if (Rank==0) cerr << "[grid=" << gid+1 << "] invalid viscosity -> model option!" << endl;
				MPI_Abort(MPI_COMM_WORLD,-1);
			}
		}
	
		if (material_input[gid].get_double("Cp").is_found) {
			Cp_model=CONSTANT;
			Cp_value=material_input[gid].get_double("Cp");
		} else if (material_input[gid].section("Cp").is_found) {
			if (material_input[gid].section("Cp").get_string("model")=="schomate") {
				Cp_model=POLY;
				Cp_poly.set("schomate",
					   material_input[gid].section("Cp").get_int("numberofpieces"),
					   material_input[gid].section("Cp").get_doubleList("coefficients"));
			} else {
				if (Rank==0) cerr << "[grid=" << gid+1 << "] invalid Cp -> model option!" << endl;
				MPI_Abort(MPI_COMM_WORLD,-1);
			}
		}
		
		if (material_input[gid].get_double("thermalconductivity").is_found) {
			lambda_model=CONSTANT;
			lambda=material_input[gid].get_double("thermalconductivity");
		} else if (material_input[gid].section("thermalconductivity").is_found) {
			if (material_input[gid].section("thermalconductivity").get_string("model")=="constantPrandtl") {
				lambda_model=PRANDTL;
				Pr=material_input[gid].section("thermalconductivity").get_double("Pr");
			} else {
				if (Rank==0) cerr << "[grid=" << gid+1 << "] invalid thermalconductivity -> model option!" << endl;
				MPI_Abort(MPI_COMM_WORLD,-1);
			}
		}

	} else {
		// If a material file is not specified
		// Ideal gas (or solid) with constant properties is assumed
		eos_model=IDEAL_GAS;
		Mw=input.section("grid",gid).subsection("material").get_double("molarmass");
		R=UNIV_GAS_CONST/Mw;
		density=input.section("grid",gid).subsection("material").get_double("density");
		gamma=input.section("grid",gid).subsection("material").get_double("gamma");
		visc_model=CONSTANT;
		mu=input.section("grid",gid).subsection("material").get_double("viscosity");
		lambda_model=CONSTANT;
		lambda=input.section("grid",gid).subsection("material").get_double("thermalconductivity");
		Cp_model=CONSTANT;
		Cp_value=input.section("grid",gid).subsection("material").get_double("Cp");
	}
	Pref=input.section("reference").get_double("p");
  	Tref=input.section("reference").get_double("T");
}

double MATERIAL::rho (double p, double T) { 
	return (p+Pref)/(R*(T+Tref));
}

double MATERIAL::p (double rho, double T) {
	return rho*R*(T+Tref)-Pref;
}

double MATERIAL::T (double p, double rho) {
	return (p+Pref)/(R*rho)-Tref;
}

double MATERIAL::a (double p, double T) { 
	return sqrt(gamma*(p+Pref)/rho(p,T));
}

double MATERIAL::viscosity (double T) {
	if (visc_model==CONSTANT) {
		return mu;
	} else if (visc_model==SUTHERLANDS) {
		return sut_mu_ref*pow((T+Tref)/sut_T_ref,1.5)*(sut_T_ref+sut_S)/(T+Tref+sut_S);
	}
}

double MATERIAL::therm_cond (double T) {
	if (lambda_model==CONSTANT) {
		return lambda;
	} else if (lambda_model==PRANDTL) {
		return Cp(T)*viscosity(T)/Pr;
	}
}

double MATERIAL::Cp (double T) {
	if (eos_model==IDEAL_GAS) {
		return gamma*R/(gamma-1.);
	} else {
		if (Cp_model==CONSTANT) {
			return Cp_value;
		} else if (Cp_model==POLY) {
			return Cp_poly.eval(T+Tref);
		}
	}
}
