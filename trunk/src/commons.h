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
#ifndef COMMONS_H
#define COMMONS_H
#include <stdlib.h>
#include <string>
#include <vector>
#include "grid.h"
#include "eos.h"

#define UNIV_GAS_CONST 8.31447215

#define NONE -1
// Options of EQUATIONS
#define NS 1
#define EULER 2
// Options for TURBULENCE_MODEL
#define KOMEGA 1
// Options for TIME_INTEGRATOR
#define FORWARD_EULER 1
#define BACKWARD_EULER 2
// Options for TIME_STEP_TYPE
#define FIXED 1
#define CFL_MAX 2
#define CFL_LOCAL 3
// Options for CONVECTIVE_FLUX_FUNCTION
#define ROE 1
#define AUSM_PLUS_UP 2
// Options for LIMITER 
#define MINMOD 1
#define DOUBLEMINMOD 2
#define HARMONIC 3
#define SUPERBEE 4
// Options for PRECONDITIONER
#define WS95 1
// Options for OUTPUT_FORMAT
#define TECPLOT 1
#define VTK 2
// Options for Boundary Condition Types
#define SYMMETRY 1
#define SLIP 2
#define NOSLIP 3
#define INLET 4
#define OUTLET 5
// Options for Boundary Condition Type variants (kind)
#define FIXED_PRESSURE 1
#define FIXED_PRESSURE_ENTRAINMENT 2
// Face bc type numbering
#define INTERNAL -1
#define UNASSIGNED -2
#define GHOST -3
// Options for order
#define FIRST 1
#define SECOND 2
// Options for thermal boundary condition
#define FIXED_T 1
#define FIXED_Q 2
#define ADIABATIC 3
// Options for Equation of State
#define IDEAL_GAS 1
// Options for BC specifications
#define BC_RHO 1
#define BC_P 2
#define BC_T 3
#define BC_STATE 4
#define BC_V 5
#define BC_EXTRAPOLATE 6
#define BC_REFLECT 7

extern int EQUATIONS,TURBULENCE_MODEL;
extern int TIME_INTEGRATOR,TIME_STEP_TYPE;
extern int CONVECTIVE_FLUX_FUNCTION,CONVECTIVE_FLUX_FUNCTION_JAC,LIMITER,PRECONDITIONER;
extern int OUTPUT_FORMAT;

extern int Rank,np;
extern int nSolVar;
extern double dt,dtTarget,CFLmax,CFLmaxTarget,CFLlocal,CFLlocalTarget;
extern int timeStep,restart;
extern double Minf,Pref,Tref;
extern int order;
extern double limiter_sharpening;
extern int jacobianUpdateFreq;
extern double Gamma,gmp1,gmm1,viscosity,conductivity;
extern double eosType,molarMass;
extern int outFreq, restartFreq;

extern double sqrt_machine_error;
extern bool ramp; 
extern double ramp_initial,ramp_growth;

extern int probeFreq,integrateBoundaryFreq;
extern int bcCount;

extern Grid grid;
extern EOS eos;

// Iterators
extern std::vector<Cell>::iterator cit;
extern std::vector<Node>::iterator nit;
extern std::vector<int>::iterator it;
extern std::vector<int>::iterator it2;
extern std::vector<double>::iterator dit;

#endif