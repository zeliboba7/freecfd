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
#include "commons.h"
#include "inputs.h"

void check_inputs(InputFile &input) {
	
	// While the input class is awesome, it is slow for frequently accessed variables
	// Pass those to global variables (defined in commons.h) and do some sanity check in the mean time
	
	string option;
	option=input.get_string("equations");
	if (option=="NS") {
		EQUATIONS=NS;
	} else if (option=="Euler") {
		EQUATIONS=EULER;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry equations=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are NS and Euler" << endl;
			exit(1);
		}
	}
		
	option=input.get_string("turbulenceModel");
	if (option=="none") {
		TURBULENCE_MODEL=NONE;
	} else if (option=="k-omega") {
		TURBULENCE_MODEL=KOMEGA;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry turbulenceModel=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are none and k-omega" << endl;
			exit(1);
		}
	}
	
	option=input.section("timeMarching").get_string("integrator");
	if (option=="forwardEuler") {
		TIME_INTEGRATOR=FORWARD_EULER;
	} else if (option=="backwardEuler") {
		TIME_INTEGRATOR=BACKWARD_EULER;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry timeMarching -> integrator=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are forwardEuler and backwardEuler" << endl;
			exit(1);
		}
	}
	
	TIME_STEP_TYPE=NONE;
	if (input.section("timeMarching").get_double("stepSize").is_found) {
		TIME_STEP_TYPE=FIXED; dt=input.section("timeMarching").get_double("stepSize");
	}
	if (input.section("timeMarching").get_double("CFLmax").is_found) {
		if (TIME_STEP_TYPE==FIXED) {
			cerr << "[E] Input entry timeMarching -> CFLmax can't be specified together with stepSize!!" << endl;
			exit(1);
		}
		TIME_STEP_TYPE=CFL_MAX; CFLmax=input.section("timeMarching").get_double("CFLmax");
	}
	else if (input.section("timeMarching").get_double("CFLlocal").is_found) {
		if (TIME_STEP_TYPE==FIXED) {
			cerr << "[E] Input entry timeMarching -> CFLlocal can't be specified together with stepSize!!" << endl;
			exit(1);
		}
		if (TIME_STEP_TYPE==CFL_MAX) {
			cerr << "[E] Input entry timeMarching -> CFLlocal can't be specified together with CFLmax!!" << endl;
			exit(1);
		}
		TIME_STEP_TYPE=CFL_LOCAL; CFLmax=input.section("timeMarching").get_double("CFLlocal");
	}
	
	option=input.section("numericalOptions").get_string("convectiveFlux");
	if (option=="ROE") {
		CONVECTIVE_FLUX_FUNCTION=ROE;
	} else if (option=="AUSM+up") {
		CONVECTIVE_FLUX_FUNCTION=AUSM_PLUS_UP;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry numericalOptions -> convectiveFlux=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are Roe and AUSM+up" << endl;
			exit(1);
		}
	}
	
	option=input.section("numericalOptions").get_string("preconditioner");
	if (option=="none") {
		PRECONDITIONER=NONE;
	} else if (option=="ws95") {
		PRECONDITIONER=WS95;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry numericalOptions -> preconditioner=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are none and ws95" << endl;
			exit(1);
		}
	}
	
	option=input.section("numericalOptions").get_string("order");
	if (option=="first") {
		order=1;
	} else if (option=="second") {
		order=2;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry numericalOptions -> order=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are first and second" << endl;
			exit(1);
		}
	}
	
	option=input.section("numericalOptions").get_string("limiter");
	if (option=="none") {
		LIMITER=NONE;
	} else if (option=="minmod") {
		LIMITER=MINMOD;
	} else if (option=="superbee") {
		LIMITER=SUPERBEE;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry numericalOptions -> limiter=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are none, minmod and superbee" << endl;
			exit(1);
		}
	}
	
	option=input.section("writeOutput").get_string("format");
	if (option=="tecplot") {
		OUTPUT_FORMAT=TECPLOT;
	} else if (option=="vtk") {
		OUTPUT_FORMAT=VTK;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry writeOutput -> format=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are tecplot and vtk" << endl;
			exit(1);
		}
	}
	
	viscosity=input.section("fluidProperties").get_double("viscosity");
	Gamma=input.section("fluidProperties").get_double("gamma");
	Minf=input.section("reference").get_double("Mach");
	Pref=input.section("reference").get_double("p");
	jacobianUpdateFreq=input.section("jacobian").get_int("updateFrequency");
	outFreq=input.section("writeOutput").get_int("plotFrequency");
	restartFreq=input.section("writeOutput").get_int("restartFrequency");
	dt=input.section("timeMarching").get_double("stepSize");
	dtTarget=dt;
	CFLmax=input.section("timeMarching").get_double("CFLmax");
	CFLmaxTarget=CFLmax;
	CFLlocal=input.section("timeMarching").get_double("CFLlocal");
	CFLlocalTarget=CFLlocal;
	ramp=input.section("timeMarching").subsection("ramp").is_found;
	ramp_initial=input.section("timeMarching").subsection("ramp").get_double("initial");
	ramp_growth=input.section("timeMarching").subsection("ramp").get_double("growth");
	probeFreq=input.section("probes").get_int("frequency");
	loadFreq=input.section("loads").get_int("frequency");
	bcCount=input.section("boundaryConditions").subsection("BC",0).count;
	
}