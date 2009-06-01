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

int DIMENSION,EQUATIONS,TURBULENCE_MODEL,TURBULENCE_FILTER;
int TIME_INTEGRATOR,TIME_STEP_TYPE;
int PS_TIME_STEP_TYPE;
int CONVECTIVE_FLUX_FUNCTION,CONVECTIVE_FLUX_FUNCTION_JAC,LIMITER,PRECONDITIONER;
int OUTPUT_FORMAT;
bool FLAMELET,GRAD_TEST;

int Rank,np;
double omegaLowLimit,kLowLimit,kHighLimit,viscosityRatioLimit,turbulenceFilterSize;
double dt_relax,dt_current,dt_target,dt_min,dt_max,CFLmax,CFLmaxTarget,CFLlocal,CFLlocalTarget;
int timeStep,restart;
int ps_timeStep;
double ps_dt_relax,ps_dt_current,ps_dt_min,ps_dt_max,ps_CFLmax,ps_CFLlocal;
double Minf,Pref,Tref;
int order;
double limiter_sharpening;
int jacobianUpdateFreq;
double Gamma,gmp1,gmm1,viscosity,conductivity;
double eosType,molarMass;
int outFreq, restartFreq;

double sqrt_machine_error;
bool ramp; 
double ramp_initial,ramp_growth;

int probeFreq,integrateBoundaryFreq;
int bcCount;

Grid grid;
EOS eos;

// Iterators
std::vector<Cell>::iterator cit;
std::vector<Cell>::iterator cit2;
std::vector<Node>::iterator nit;
std::vector<int>::iterator it;
std::vector<int>::iterator it2;
std::vector<double>::iterator dit;

// TODO move these over to maps class
std::vector< std::vector<int> > sendCells;
std::vector<int>  recvCount;
