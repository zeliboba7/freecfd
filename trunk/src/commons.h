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
#include <string>
#include <vector>
#include "grid.h"

#define NONE -1
// Options of EQUATIONS
#define NS 0
#define EULER 1
// Options for TURBULENCE_MODEL
#define KOMEGA 0
// Options for TIME_INTEGRATOR
#define FORWARD_EULER 0
#define BACKWARD_EULER 1
// Options for TIME_STEP_TYPE
#define FIXED 0
#define CFL_MAX 1
#define CFL_LOCAL 2
// Options for CONVECTIVE_FLUX_FUNCTION
#define ROE 1
#define AUSM_PLUS_UP 1
// Options for LIMITER 
#define MINMOD 0
#define SUPERBEE 1
// Options for PRECONDITIONER
#define WS95 0
// Options for OUTPUT_FORMAT
#define TECPLOT 0
#define VTK 1
// Options for Boundary Condition Types
#define SYMMETRY 0
#define SLIP 1
#define NOSLIP 2
#define INLET 3
#define OUTLET 4
// Options for Boundary Condition Type variants (kind)
#define FIXED_PRESSURE 0
#define FIXED_PRESSURE_ENTRAINMENT 1

extern int EQUATIONS,TURBULENCE_MODEL;
extern int TIME_INTEGRATOR,TIME_STEP_TYPE;
extern int CONVECTIVE_FLUX_FUNCTION,LIMITER,PRECONDITIONER;
extern int OUTPUT_FORMAT;

extern int Rank,np;
extern double dt,dtTarget,CFLmax,CFLmaxTarget,CFLlocal,CFLlocalTarget;
extern int timeStep,restart;
extern double Minf, Pref;
extern int order;
extern double limiter_sharpening;
extern int jacobianUpdateFreq;
extern double Gamma,viscosity;
extern int outFreq, restartFreq;

extern double sqrt_machine_error;
extern bool ramp; 
extern double ramp_initial,ramp_growth;

extern int probeFreq,loadFreq;
extern int bcCount;

extern Grid grid;
// Iterators
extern std::vector<Cell>::iterator cit;

#endif