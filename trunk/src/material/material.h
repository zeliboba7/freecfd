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
#ifndef MATERIAL_H
#define MATERIAL_H

#define UNIV_GAS_CONST 8314.47215
// Options for Equation of State
#define IDEAL_GAS 1
// Options for viscosity
#define CONSTANT 1
#define SUTHERLANDS 2
// Options for thermal conductivity
#define PRANDTL 2
// Options for Cp
#define POLY 2

#include <cmath>
#include "inputs.h"
#include "polynomial.h"
#include "commons.h"
extern InputFile input;
extern vector<InputFile> material_input;

class MATERIAL {
public:
	int eos_model;
	int visc_model,lambda_model,Cp_model;
	double Mw; // Molar mass
	double R; // Ideal gas constant
	double gamma;
	double Cp_value;
	Polynomial Cp_poly;
	double density;
	double Pref,Tref; // Reference properties (everything is relative to these)
	// Constant viscosity and thermal conductivity variables
	double mu, lambda;
	// Sutherlands model variables
	double sut_mu_ref,sut_T_ref,sut_S;
	// Prandtl number
	double Pr;
	
	MATERIAL(); // Constructor
	void set(int gid);
	// Equation of State functions
	double rho (double p, double T);
	double p (double rho, double T);
	double T (double p, double rho);
	double a (double p, double T);
	// Viscosity 
	double viscosity (double T=0.);
	double therm_cond (double T=0.);
	double Cp (double T=0.);
};

#endif
