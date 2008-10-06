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
#include "inputs.h"

void read_inputs(InputFile &input) {

	input.register_section("equations");
	input.section["equations"].register_string("set");
	// defaults
	input.section["equations"].strings["set"]="Euler";
	input.read_section("equations");

	input.register_section("timeMarching");
	input.section["timeMarching"].register_string("integrator");
	input.section["timeMarching"].register_string("preconditioner");
	input.section["timeMarching"].register_string("type");
	input.section["timeMarching"].register_double("step");
	input.section["timeMarching"].register_double("CFL");
	input.section["timeMarching"].register_int("numberOfSteps");
	// defaults
	input.section["timeMarching"].strings["type"]="CFL";
	input.section["timeMarching"].doubles["CFL"]=1.;
	input.section["timeMarching"].ints["numberOfSteps"]=100;
	input.read_section("timeMarching");

	input.register_section("initialConditions");
	input.section["initialConditions"].register_numberedSubsection("region");
	input.section["initialConditions"].numberedSubsections["region"].register_string("type");
	input.section["initialConditions"].numberedSubsections["region"].register_Vec3D("box_1");
	input.section["initialConditions"].numberedSubsections["region"].register_Vec3D("box_2");
	input.section["initialConditions"].numberedSubsections["region"].register_double("rho");
	input.section["initialConditions"].numberedSubsections["region"].register_Vec3D("center");
	input.section["initialConditions"].numberedSubsections["region"].register_double("radius");
	input.section["initialConditions"].numberedSubsections["region"].register_Vec3D("v");
	input.section["initialConditions"].numberedSubsections["region"].register_double("p");
	input.section["initialConditions"].numberedSubsections["region"].register_double("k");
	input.section["initialConditions"].numberedSubsections["region"].register_double("omega");
	input.read_section("initialConditions");

	input.register_section("boundaryConditions");
	input.section["boundaryConditions"].register_numberedSubsection("BC");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_string("type");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_string("kind");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_string("region");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_Vec3D("box_1");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_Vec3D("box_2");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_string("pick");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_double("rho");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_double("p");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_double("k");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_double("omega");
	input.section["boundaryConditions"].numberedSubsections["BC"].register_Vec3D("v");
	input.read_section("boundaryConditions");

	input.register_section("fluidProperties");
	input.section["fluidProperties"].register_double("gamma");
	input.section["fluidProperties"].register_double("Pref");
	input.section["fluidProperties"].register_subsection("viscosity");
	input.section["fluidProperties"].subsections["viscosity"].register_string("type");
	input.section["fluidProperties"].subsections["viscosity"].register_double("value");
	// defaults
	input.section["fluidProperties"].subsections["viscosity"].strings["type"]="fixed";
	input.read_section("fluidProperties");

	input.register_section("numericalOptions");
	input.section["numericalOptions"].register_string("flux");
	input.section["numericalOptions"].register_string("order");
	input.section["numericalOptions"].register_string("limiter");
	input.section["numericalOptions"].register_double("Minf");
	input.section["numericalOptions"].register_double("sharpeningFactor");
	// defaults
	input.section["numericalOptions"].strings["order"]="first";
	input.section["numericalOptions"].strings["limiter"]="superbee";
	input.section["numericalOptions"].doubles["sharpeningFactor"]=0.;
	input.read_section("numericalOptions");

	input.register_section("linearSolver");
	input.section["linearSolver"].register_double("relTolerance");
	input.section["linearSolver"].register_double("absTolerance");
	input.section["linearSolver"].register_int("maxIterations");
	input.read_section("linearSolver");
	
	input.register_section("output");
	input.section["output"].register_string("format");
	input.section["output"].register_int("outFreq");
	input.section["output"].register_int("restartFreq");
	input.read_section("output");

	input.register_section("probes");
	input.section["probes"].register_int("frequency");
	input.section["probes"].register_numberedSubsection("probe");
	input.section["probes"].numberedSubsections["probe"].register_Vec3D("coord");
	input.read_section("probes");

	input.register_section("loads");
	input.section["loads"].register_int("frequency");
	input.section["loads"].register_numberedSubsection("load");
	input.section["loads"].numberedSubsections["load"].register_int("bc");
	input.read_section("loads");

	input.register_section("turbulence");
	input.section["turbulence"].register_string("model");
	input.read_section("turbulence");

	
}
