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
#include "grid.h"
#include "inputs.h"
extern bool grad_test;// DEBUG

bool within_box(Vec3D centroid, Vec3D box_1, Vec3D box_2);

void initialize(Grid &grid, InputFile &input) {

	// Loop through each initial condition region and apply sequentially
	numberedSubsection region=input.section["initialConditions"].numberedSubsections["region"];
	for (unsigned int r=0;r<region.count;++r) {
		if (region.strings[r]["type"]=="box") {
			for (unsigned int c = 0;c < grid.cellCount;++c) {
				Vec3D box_1=region.Vec3Ds[r]["box_1"];
				Vec3D box_2=region.Vec3Ds[r]["box_2"];
				if (within_box(grid.cell[c].centroid,box_1,box_2)) {
					// The cell centroid is inside the box region
					grid.cell[c].rho = region.doubles[r]["rho"];
					grid.cell[c].v = region.Vec3Ds[r]["v"];
					grid.cell[c].p = region.doubles[r]["p"];
					grid.cell[c].k = region.doubles[r]["k"];
					grid.cell[c].omega = region.doubles[r]["omega"];
				}
			}
		} else if (region.strings[r]["type"]=="circle") {
			double radius=region.doubles[r]["radius"];
			Vec3D center=region.Vec3Ds[r]["center"];
			Vec3D zComp=center;
			Vec3D center2cell;
			zComp.comp[0]=0.;zComp.comp[1]=0.;
			center=center-zComp;
			for (unsigned int c = 0;c < grid.cellCount;++c) {
				zComp=grid.cell[c].centroid;  zComp.comp[0]=0.;zComp.comp[1]=0.;
				center2cell= (grid.cell[c].centroid-zComp)-center;
				if (fabs(center2cell) <=radius) {
					// The cell centroid is inside the box region
					grid.cell[c].rho = region.doubles[r]["rho"];
					grid.cell[c].v = region.Vec3Ds[r]["v"];
					grid.cell[c].v.comp[0]= (region.Vec3Ds[r]["v"].comp[0]*center2cell).comp[0];
					grid.cell[c].v.comp[1]= (region.Vec3Ds[r]["v"].comp[0]*center2cell).comp[1];
					grid.cell[c].p = region.doubles[r]["p"];
					grid.cell[c].k = region.doubles[r]["k"];
					grid.cell[c].omega = region.doubles[r]["omega"];
				}
			}
		} else if (region.strings[r]["type"]=="sphere") {
			double radius=region.doubles[r]["radius"];
			Vec3D center=region.Vec3Ds[r]["center"];
			Vec3D center2cell,center2cellUnit;
			Vec3D unitX, unitY, unitZ;
			unitX.comp[0]=1.; unitX.comp[1]=0.;unitX.comp[2]=0.;
			unitY.comp[1]=0.; unitY.comp[1]=1.;unitY.comp[2]=0.;
			unitZ.comp[2]=0.; unitZ.comp[1]=0.;unitZ.comp[2]=1.;
			for (unsigned int c = 0;c < grid.cellCount;++c) {
				center2cell= grid.cell[c].centroid-center;
				center2cellUnit=center2cell/fabs(center2cell);
				if (fabs(center2cell) <=radius) {
					// The cell centroid is inside the box region
					grid.cell[c].rho = region.doubles[r]["rho"];
					grid.cell[c].v = region.Vec3Ds[r]["v"].comp[0]*center2cellUnit;
					grid.cell[c].p = region.doubles[r]["p"];
					grid.cell[c].k = region.doubles[r]["k"];
					grid.cell[c].omega = region.doubles[r]["omega"];
				}
			}
		}
	}

	bool mu_is_fixed=false;
	if (input.section["fluidProperties"].subsections["viscosity"].strings["type"]=="fixed") mu_is_fixed=true;
	double mu_fixed=input.section["fluidProperties"].subsections["viscosity"].doubles["value"];
	for (unsigned int c=0;c<grid.cellCount;++c) {
		if (mu_is_fixed) {
			grid.cell[c].mu=mu_fixed;
		}
		for (unsigned int i=0;i<7;++i) {
			grid.cell[c].flux[i]=0.;
		}
	}

	if (grad_test) { // DEBUG
		for (unsigned int c=0;c<grid.cellCount;++c) { // DEBUG
			//grid.cell[c].rho=2.*grid.cell[c].centroid.comp[0]+2.; // DEBUG
			grid.cell[c].rho=2.*grid.cell[c].centroid.comp[0]*grid.cell[c].centroid.comp[0]+2.; // DEBUG
			//grid.cell[c].rho=1.; // DEBUG
		}// DEBUG
	}// DEBUG
	
	return;
}
