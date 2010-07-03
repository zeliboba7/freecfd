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
#include "grid.h"
#include "inputs.h"
#include "variable.h"
#include "ns.h"

extern InputFile input;
extern vector<Grid> grid;
extern vector<NavierStokes> ns;
// Time step for each grid
extern vector<Variable<double> > dt;

#define NONE -1
// Options for time_integrator
#define FORWARD_EULER 1
#define BACKWARD_EULER 2
// Options for time_step_type
#define FIXED 1
#define CFL_MAX 2
#define CFL_LOCAL 3

int time_integrator, time_step_type;
double time_step_current,time_step_target;
double CFLmax,CFLmaxTarget,CFLlocal,CFLlocalTarget;
bool time_step_ramp;
double time_step_ramp_initial,time_step_ramp_growth;
