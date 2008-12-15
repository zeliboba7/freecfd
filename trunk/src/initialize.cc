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
#include "grid.h"
#include "inputs.h"
extern bool grad_test;// DEBUG

bool withinBox(Vec3D point, Vec3D corner_1, Vec3D corner_2);
bool withinCylinder(Vec3D point, Vec3D center, double radius, Vec3D axisDirection, double height);
bool withinSphere(Vec3D point, Vec3D center, double radius);

void initialize(InputFile &input) {

	// Loop through each initial condition region and apply sequentially
	int count=input.section("initialConditions").subsection("IC",0).count;
	for (unsigned int ic=0;ic<count;++ic) {
		// Store the reference to current IC region
		Subsection &region=input.section("initialConditions").subsection("IC",ic);
		Vec3D regionV=region.get_Vec3D("v");
		// If region is specified with a box method
		if (region.get_string("region")=="box") {
			// Loop the cells
			for (unsigned int c=0;c<grid.cellCount;++c) {
				// Check if the cell centroid is inside the box region
				if (withinBox(grid.cell[c].centroid,region.get_Vec3D("corner_1"),region.get_Vec3D("corner_2"))) {
					// Assign specified values
					grid.cell[c].rho=region.get_double("rho");
					grid.cell[c].v=region.get_Vec3D("v");
					grid.cell[c].p=region.get_double("p");
					grid.cell[c].k=region.get_double("k");
					grid.cell[c].omega=region.get_double("omega");
				}
			}
		} else if (region.get_string("region")=="cylinder") {
			// Loop the cells
			for (unsigned int c=0;c<grid.cellCount;++c) {
				// Check if the cell centroid is inside the cylinder region
				Vec3D axisDirection=region.get_Vec3D("axisDirection");
				axisDirection=axisDirection.norm();
				if (withinCylinder(grid.cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"),axisDirection,region.get_double("height"))) {
					grid.cell[c].rho=region.get_double("rho");
					grid.cell[c].p=region.get_double("p");
					grid.cell[c].k=region.get_double("k");
					grid.cell[c].omega=region.get_double("omega");
					// first component of the specified velocity is interpreted as the axial velocity
					// second component of the specified velocity is interpreted as the radial velocity
					// third component of the specified velocity is interpreted as the circumferential velocity
					Vec3D radialPoint=grid.cell[c].centroid-region.get_Vec3D("center");
					radialPoint=radialPoint.dot(axisDirection);	
					grid.cell[c].v=		
						regionV[0]*axisDirection /* axial component */
						+regionV[1]*radialPoint.norm() /* circumferential component */
						+regionV[2]*radialPoint.cross(axisDirection); /* circumferential component */
				}
			}
		} else if (region.get_string("region")=="sphere") {
			// Loop the cells
			for (unsigned int c=0;c<grid.cellCount;++c) {
				// Check if the cell centroid is inside the sphere region
				if (withinSphere(grid.cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"))) {
					grid.cell[c].rho=region.get_double("rho");
					grid.cell[c].p=region.get_double("p");
					grid.cell[c].k=region.get_double("k");
					grid.cell[c].omega=region.get_double("omega");
					// first component of the specified velocity is interpreted as the radial velocity
					// second and third components of the specified velocity are ignored
					grid.cell[c].v=	regionV[0]*(grid.cell[c].centroid-region.get_Vec3D("center"));
				}
			}
		}
	}

	for (unsigned int c=0;c<grid.cellCount;++c) {
		grid.cell[c].mu=viscosity;
		for (unsigned int i=0;i<7;++i) {
			grid.cell[c].flux[i]=0.;
			grid.cell[c].update[i]=0.;
		}
	}
	
	if (grad_test) { // DEBUG
		for (unsigned int c=0;c<grid.cellCount;++c) { // DEBUG
			//grid.cell[c].rho=2.*grid.cell[c].centroid.comp[0]+2.; // DEBUG
			grid.cell[c].rho=2.*grid.cell[c].centroid.comp[0]*grid.cell[c].centroid.comp[0]+2.*grid.cell[c].centroid.comp[0]+2.; // DEBUG
			//grid.cell[c].rho=1.; // DEBUG
		}// DEBUG
	}// DEBUG
	
	
	
	return;
}
