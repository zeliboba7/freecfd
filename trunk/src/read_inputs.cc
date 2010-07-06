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
#include "inputs.h"

extern InputFile input;
extern vector<InputFile> material_input;

void read_inputs(void) {

	// These will make more sense as you read the input registration lines below
	bool required=true; bool optional=false;
	bool numbered=true; bool single=false;
	
	input.registerSection("reference",single,optional);
	input.section("reference").register_double("p",optional,0.);
	input.section("reference").register_double("T",optional,0.);
	input.section("reference").register_double("Mach",optional,1.);
	input.read("reference");
	
	input.registerSection("grid",numbered,required);
	input.section("grid",0).register_string("file",required);
	input.section("grid",0).register_int("dimension",optional,3);
	
	input.section("grid",0).register_string("equations",required);
	
	input.section("grid",0).registerSubsection("writeoutput",single,required);
	input.section("grid",0).subsection("writeoutput").register_int("plotfrequency",optional,1e20);
	input.section("grid",0).subsection("writeoutput").register_int("restartfrequency",optional,1e20);
	input.section("grid",0).subsection("writeoutput").register_stringList("variables",required);

	input.section("grid",0).register_string("material",optional,"none");
	input.section("grid",0).registerSubsection("material",single,optional);
	// Default air values assigned
	input.section("grid",0).subsection("material").register_double("gamma",optional,1.4);
	input.section("grid",0).subsection("material").register_double("molarmass",optional,28.97);
	input.section("grid",0).subsection("material").register_double("viscosity",optional,0.);
	input.section("grid",0).subsection("material").register_double("thermalconductivity",optional,0.);
	input.section("grid",0).subsection("material").register_double("Cp",optional,1000.);
	input.section("grid",0).subsection("material").register_double("density",optional,0.);
	
	input.section("grid",0).registerSubsection("IC",numbered,required);
	input.section("grid",0).subsection("IC",0).register_string("region",optional,"box");
	input.section("grid",0).subsection("IC",0).register_Vec3D("corner_1",optional,-1.e20);
	input.section("grid",0).subsection("IC",0).register_Vec3D("corner_2",optional,1e20);
	input.section("grid",0).subsection("IC",0).register_Vec3D("center",optional);
	input.section("grid",0).subsection("IC",0).register_Vec3D("axisdirection",optional);
	input.section("grid",0).subsection("IC",0).register_double("radius",optional);
	input.section("grid",0).subsection("IC",0).register_double("height",optional);
	input.section("grid",0).subsection("IC",0).register_double("p",optional);
	input.section("grid",0).subsection("IC",0).register_Vec3D("V",optional);
	input.section("grid",0).subsection("IC",0).register_double("T",optional);
	input.section("grid",0).subsection("IC",0).register_double("rho",optional);
	input.section("grid",0).subsection("IC",0).register_double("turbulenceintensity",optional,1.e-2);
	input.section("grid",0).subsection("IC",0).register_double("eddyviscosityratio",optional,0.1);
	
	input.section("grid",0).registerSubsection("BC",numbered,required);
	input.section("grid",0).subsection("BC",0).register_string("type",required);
	input.section("grid",0).subsection("BC",0).register_string("kind",optional,"none");
	input.section("grid",0).subsection("BC",0).register_string("region",optional,"gridfile");
	input.section("grid",0).subsection("BC",0).register_string("interface",optional,"none");	
	input.section("grid",0).subsection("BC",0).register_Vec3D("corner_1",optional);
	input.section("grid",0).subsection("BC",0).register_Vec3D("corner_2",optional);
	input.section("grid",0).subsection("BC",0).register_string("pick",optional,"override");
	input.section("grid",0).subsection("BC",0).register_double("p",optional);
	input.section("grid",0).subsection("BC",0).register_double("mdot",optional);
	input.section("grid",0).subsection("BC",0).register_double("qdot",optional);
	input.section("grid",0).subsection("BC",0).register_Vec3D("V",optional);
	input.section("grid",0).subsection("BC",0).register_double("T",optional);
	input.section("grid",0).subsection("BC",0).register_double("rho",optional);
	input.section("grid",0).subsection("BC",0).register_double("k",optional,0.);
	input.section("grid",0).subsection("BC",0).register_double("omega",optional,0.);
	
	input.section("grid",0).registerSubsection("turbulence",single,optional);
	input.section("grid",0).subsection("turbulence").register_double("relativetolerance",optional,1.e-6);
	input.section("grid",0).subsection("turbulence").register_double("absolutetolerance",optional,1.e-12);
	input.section("grid",0).subsection("turbulence").register_int("maximumiterations",optional,10);	
	input.section("grid",0).subsection("turbulence").register_string("model",optional,"none");
	input.section("grid",0).subsection("turbulence").register_string("order",optional,"second");
	
	input.section("grid",0).registerSubsection("navierstokes",single,optional);
	input.section("grid",0).subsection("navierstokes").register_double("relativetolerance",optional,1.e-6);
	input.section("grid",0).subsection("navierstokes").register_double("absolutetolerance",optional,1.e-12);
	input.section("grid",0).subsection("navierstokes").register_int("maximumiterations",optional,10);	
	input.section("grid",0).subsection("navierstokes").register_string("limiter",optional,"none");
	input.section("grid",0).subsection("navierstokes").register_double("limiter_threshold",optional,1.);
	input.section("grid",0).subsection("navierstokes").register_string("order",optional,"second");
	input.section("grid",0).subsection("navierstokes").register_string("convectiveflux",optional,"AUSM+up");
	
	input.section("grid",0).registerSubsection("transform",numbered,optional);
	input.section("grid",0).subsection("transform",0).register_string("function",optional);
	input.section("grid",0).subsection("transform",0).register_Vec3D("anchor",optional,0.);
	input.section("grid",0).subsection("transform",0).register_Vec3D("begin",optional,0.);
	input.section("grid",0).subsection("transform",0).register_Vec3D("end",optional,0.);
	input.section("grid",0).subsection("transform",0).register_double("factor",optional,0.);
	input.section("grid",0).subsection("transform",0).register_Vec3D("axis",optional,0.);
	input.section("grid",0).subsection("transform",0).register_double("angle",optional,0.);
	
	input.read("grid",0);
	
	input.registerSection("timemarching",single,required);
	input.section("timemarching").register_string("integrator",optional,"backwardEuler");
	input.section("timemarching").register_double("stepsize",optional,1.);
	input.section("timemarching").register_double("CFLmax",optional,1000.);
	input.section("timemarching").register_double("CFLlocal",optional,1000.);
	input.section("timemarching").registerSubsection("ramp",single,optional);
	input.section("timemarching").subsection("ramp").register_double("initial",optional,1.);
	input.section("timemarching").subsection("ramp").register_double("growth",optional,1.2);
	input.section("timemarching").register_int("numberofsteps",required);
	input.read("timemarching");

	input.readEntries();
	
	// Read the material file for each grid
	material_input.resize(input.section("grid",0).count);
	for (int gid=0;gid<material_input.size();++gid) {
		// If a material specification exists
		if(input.section("grid",0).get_string("material").is_found) {
			// Default values are set to air properties
			string fileName=input.section("grid",0).get_string("material");
			fileName+=".mat";
			material_input[gid].setFile(fileName);
			material_input[gid].register_double("molar mass",optional,28.97);
			material_input[gid].register_double("density",optional,0.);
			material_input[gid].register_double("gamma",optional,1.4);
			material_input[gid].registerSection("equationofstate",single,optional);
			material_input[gid].section("equationofstate").register_string("model",optional,"idealgas");
			material_input[gid].read("equationofstate");
			material_input[gid].register_double("viscosity",optional,1.716e-5);
			material_input[gid].registerSection("viscosity",single,optional);
			material_input[gid].section("viscosity").register_string("model",optional,"sutherlands");
			material_input[gid].section("viscosity").register_double("referenceviscosity",optional,1.716e-5);
			material_input[gid].section("viscosity").register_double("referencetemperature",optional,273.15);
			material_input[gid].section("viscosity").register_double("sutherlandtemperature",optional,110.4);			
			material_input[gid].read("viscosity");
			material_input[gid].register_double("Cp",optional,1000.);
			material_input[gid].registerSection("Cp",single,optional);
			material_input[gid].section("Cp").register_int("numberofpieces",optional);
			material_input[gid].section("Cp").register_string("model",optional,"shomate");
			material_input[gid].section("Cp").register_doubleList("coefficients",optional);
			material_input[gid].read("Cp");
			material_input[gid].register_double("thermalconductivity",optional,0.024);
			material_input[gid].registerSection("thermalconductivity",single,optional);
			material_input[gid].section("thermalconductivity").register_string("model",optional,"constantPrandtl");
			material_input[gid].section("thermalconductivity").register_double("Pr",optional,0.7);
			material_input[gid].read("thermalconductivity");
			material_input[gid].readEntries();
		}
	}

	return;
}
