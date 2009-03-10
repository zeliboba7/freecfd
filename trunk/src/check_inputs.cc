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
#include "inputs.h"

void check_inputs(InputFile &input) {
	
	// While the input class is awesome, it is slow for frequently accessed variables
	// Pass those to global variables (defined in commons.h) and do some sanity check in the mean time

	if (input.section("grid").get_int("dimension").is_found) {
		if (input.section("grid").get_int("dimension")>3 || input.section("grid").get_int("dimension")<1) {
			cerr << "[E] Input entry grid -> dimension is not recognized" << endl;
			cerr << "[E] Acceptable options are 1, 2 or 3" << endl;
			exit(1);
		}
		DIMENSION=input.section("grid").get_int("dimension");
	}

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
	
	option=input.section("turbulence").get_string("model");
	
	if (option=="none") {
		TURBULENCE_MODEL=NONE;
	} else if (option=="k-omega") {
		TURBULENCE_MODEL=KOMEGA;
	} else if (option=="k-epsilon") {
		TURBULENCE_MODEL=KEPSILON;
	} else if (option=="bsl") {
		TURBULENCE_MODEL=BSL;
	} else if (option=="sst") {
		TURBULENCE_MODEL=SST;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry turbulence -> model=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are none, k-omega, k-epsilon, bsl and sst" << endl;
			exit(1);
		}
	}
	TURBULENCE_FILTER=NONE;
	if (input.section("turbulence").subsection("filter").is_found) {
		option=input.section("turbulence").subsection("filter").get_string("type");
		if (option=="uniform") {
			TURBULENCE_FILTER=UNIFORM;
			turbulenceFilterSize=input.section("turbulence").subsection("filter").get_double("size");
		} else if(option=="local") {
			TURBULENCE_FILTER=LOCAL;
			turbulenceFilterSize=input.section("turbulence").subsection("filter").get_double("size");	
		} else {
			cerr << "[E] Input entry turbulence -> filter -> kind=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are uniform and local" << endl;
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
	dt_current=input.section("timeMarching").get_double("stepSize");
	if (input.section("timeMarching").get_double("stepSize").is_found) {
		TIME_STEP_TYPE=FIXED; 
	}
	if (input.section("timeMarching").get_double("CFLmax").is_found) {
		if (TIME_STEP_TYPE==FIXED) {
			cerr << "[E] Input entry timeMarching -> CFLmax can't be specified together with stepSize!!" << endl;
			exit(1);
		}
		TIME_STEP_TYPE=CFL_MAX; CFLmax=input.section("timeMarching").get_double("CFLmax");
		CFLmaxTarget=CFLmax;
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
		TIME_STEP_TYPE=CFL_LOCAL; CFLlocal=input.section("timeMarching").get_double("CFLlocal");
		CFLlocalTarget=CFLlocal;
	}
	else if (input.section("timeMarching").get_double("adaptive").is_found) {

		if (TIME_STEP_TYPE==CFL_MAX) {
			cerr << "[E] Input entry timeMarching -> adaptive can't be specified together with CFLmax!!" << endl;
			exit(1);
		}
		if (TIME_STEP_TYPE==CFL_LOCAL) {
			cerr << "[E] Input entry timeMarching -> adaptive can't be specified together with CFLlocal!!" << endl;
			exit(1);
		}
		TIME_STEP_TYPE=ADAPTIVE; 
		dt_relax=input.section("timeMarching").get_double("adaptive");
		dt_min=input.section("timeMarching").get_double("stepSizeMin");
		dt_max=input.section("timeMarching").get_double("stepSizeMax");
	}
	
	option=input.section("numericalOptions").get_string("convectiveFlux");
	if (option=="Roe") {
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
	
	option=input.section("numericalOptions").get_string("convectiveFluxJac");
	if (option=="Roe") {
		CONVECTIVE_FLUX_FUNCTION_JAC=ROE;
	} else if (option=="AUSM+up") {
		CONVECTIVE_FLUX_FUNCTION_JAC=AUSM_PLUS_UP;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry numericalOptions -> convectiveFluxJac=" << option << " is not recognized!!" << endl;
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
	} else if (option=="doubleMinmod") {
		LIMITER=DOUBLEMINMOD;
	} else if (option=="harmonic") {
		LIMITER=HARMONIC;
	} else if (option=="superbee") {
		LIMITER=SUPERBEE;
	} else {
		if (Rank==0) {
			cerr << "[E] Input entry numericalOptions -> limiter=" << option << " is not recognized!!" << endl;
			cerr << "[E] Acceptable options are none, minmod, doubleMinmod, harmonic and superbee" << endl;
			exit(1);
		}
	}
	
	
	option=input.section("fluidProperties").get_string("eos");
	if (option!="idealGas") {
		if (Rank==0) {
			cerr << "[E] Input entry fluidProperties -> eos=" << option << " is not recognized!!" << endl;
			cerr << "[E] Currently only available option is idealGas" << endl;
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
	
	FLAMELET=input.section("flamelet").get_string("tableFile").is_found;
	omegaLowLimit=input.section("turbulence").get_double("omegaLowLimit");
	kLowLimit=input.section("turbulence").get_double("kLowLimit");
	kHighLimit=input.section("turbulence").get_double("kHighLimit");
	viscosityRatioLimit=input.section("turbulence").get_double("viscosityRatioLimit");
	viscosity=input.section("fluidProperties").get_double("viscosity");
	conductivity=input.section("fluidProperties").get_double("thermalConductivity");
	Gamma=input.section("fluidProperties").get_double("gamma");
	gmp1=Gamma+1.; gmm1=Gamma-1.;
	Minf=input.section("reference").get_double("Mach");
	Pref=input.section("reference").get_double("p");
	Tref=input.section("reference").get_double("T");
	jacobianUpdateFreq=input.section("jacobian").get_int("updateFrequency");
	outFreq=input.section("writeOutput").get_int("plotFrequency");
	restartFreq=input.section("writeOutput").get_int("restartFrequency");
	ramp=input.section("timeMarching").subsection("ramp").is_found;
	ramp_initial=input.section("timeMarching").subsection("ramp").get_double("initial");
	ramp_growth=input.section("timeMarching").subsection("ramp").get_double("growth");
	limiter_sharpening=input.section("numericalOptions").get_double("sharpeningFactor");
	probeFreq=input.section("probes").get_int("frequency");
	integrateBoundaryFreq=input.section("integrateBoundary").get_int("frequency");
	bcCount=input.section("boundaryConditions").subsection("BC",0).count;
	
	
}
