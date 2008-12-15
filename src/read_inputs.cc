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

void check_inputs(InputFile &input);

void read_inputs(InputFile &input) {

	// These will make more sense as you read the input registration lines below
	bool required=true; bool optional=false;
	bool numbered=true; bool single=false;
	
	input.register_string("equations",optional,"NS");
	input.register_string("turbulenceModel",optional,"none");
	input.readEntries();
	
	input.registerSection("grid",optional);
	input.section("grid").register_double("scaleBy",optional,1.);
	input.readSection("grid");
	
	input.registerSection("reference",optional);
	input.section("reference").register_double("Mach",optional,1.);
	input.section("reference").register_double("p",optional,0.);
	input.readSection("reference");
	
	input.registerSection("timeMarching",required);
	input.section("timeMarching").register_string("integrator",optional,"backwardEuler");
	input.section("timeMarching").register_double("stepSize",optional,1000);
	input.section("timeMarching").register_double("CFLmax",optional,1000);
	input.section("timeMarching").register_double("CFLlocal",optional,1000);
	input.section("timeMarching").registerSubsection("ramp",single,optional);
	input.section("timeMarching").subsection("ramp").register_double("initial",optional,1.);
	input.section("timeMarching").subsection("ramp").register_double("growth",optional,1.2);
	input.section("timeMarching").register_int("numberOfSteps",required);
	input.readSection("timeMarching");
	
	input.registerSection("numericalOptions",optional);
	input.section("numericalOptions").register_string("convectiveFlux",optional,"AUSM+up");
	input.section("numericalOptions").register_string("preconditioner",optional,"none");
	input.section("numericalOptions").register_string("order",optional,"second");
	input.section("numericalOptions").register_string("limiter",optional,"none");
	input.section("numericalOptions").register_double("sharpeningFactor",optional,0.);
	input.readSection("numericalOptions");
	
	input.registerSection("linearSolver",required);
	input.section("linearSolver").register_double("relTolerance",optional,1.e-6);
	input.section("linearSolver").register_double("absTolerance",optional,1.e-12);
	input.section("linearSolver").register_int("maxIterations",required);
	input.readSection("linearSolver");
	
	input.registerSection("jacobian",optional);
	input.section("jacobian").register_int("updateFrequency",optional,1);
	input.readSection("jacobian");
	
	input.registerSection("fluidProperties",required);
	input.section("fluidProperties").register_double("gamma",required);
	bool viscRequired= (input.get_string("equations")=="Euler")? false : true;
	input.section("fluidProperties").register_double("viscosity",viscRequired,0.); 
	input.readSection("fluidProperties");
	
	input.registerSection("writeOutput",required);
	input.section("writeOutput").register_string("format",optional,"tecplot");
	input.section("writeOutput").register_int("plotFrequency",required);
	input.section("writeOutput").register_int("restartFrequency",required);
	input.readSection("writeOutput");
	
	input.registerSection("probes",optional);
	input.section("probes").register_int("frequency",optional);
	input.section("probes").registerSubsection("probe",numbered,optional);
	input.section("probes").subsection("probe",0).register_Vec3D("coord",optional);
	input.readSection("probes");

	input.registerSection("loads",optional);
	input.section("loads").register_int("frequency",optional);
	input.section("loads").registerSubsection("load",numbered,optional);
	input.section("loads").subsection("load",0).register_int("bc",optional);
	input.readSection("loads");
	
	input.registerSection("initialConditions",required);
	input.section("initialConditions").registerSubsection("IC",numbered,required);
	input.section("initialConditions").subsection("IC",0).register_string("region",optional,"box");
	input.section("initialConditions").subsection("IC",0).register_Vec3D("corner_1",optional,-1.e20);
	input.section("initialConditions").subsection("IC",0).register_Vec3D("corner_2",optional,1e20);
	input.section("initialConditions").subsection("IC",0).register_Vec3D("center",optional);
	input.section("initialConditions").subsection("IC",0).register_Vec3D("axisDirection",optional);
	input.section("initialConditions").subsection("IC",0).register_double("radius",optional);
	input.section("initialConditions").subsection("IC",0).register_double("height",optional);
	input.section("initialConditions").subsection("IC",0).register_double("rho",required);
	input.section("initialConditions").subsection("IC",0).register_Vec3D("v",required);
	input.section("initialConditions").subsection("IC",0).register_double("p",required);
	input.section("initialConditions").subsection("IC",0).register_double("k",optional);
	input.section("initialConditions").subsection("IC",0).register_double("omega",optional);
	input.readSection("initialConditions");
	
	input.registerSection("boundaryConditions",required);
	input.section("boundaryConditions").registerSubsection("BC",numbered,required);
	input.section("boundaryConditions").subsection("BC",0).register_string("type",required);
	input.section("boundaryConditions").subsection("BC",0).register_string("kind",optional,"none");
	input.section("boundaryConditions").subsection("BC",0).register_string("region",optional,"gridFile");
	input.section("boundaryConditions").subsection("BC",0).register_Vec3D("corner_1",optional);
	input.section("boundaryConditions").subsection("BC",0).register_Vec3D("corner_2",optional);
	input.section("boundaryConditions").subsection("BC",0).register_string("pick",optional,"overRide");
	input.section("boundaryConditions").subsection("BC",0).register_double("rho",optional);
	input.section("boundaryConditions").subsection("BC",0).register_Vec3D("v",optional);	
	input.section("boundaryConditions").subsection("BC",0).register_double("p",optional);
	input.section("boundaryConditions").subsection("BC",0).register_double("k",optional);
	input.section("boundaryConditions").subsection("BC",0).register_double("omega",optional);
	input.readSection("boundaryConditions");
	
	check_inputs(input);
	
}