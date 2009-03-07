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
#include "grid.h"
#include "inputs.h"
#include "rans.h"
#include "flamelet.h"
extern RANS rans;
extern Flamelet flamelet;

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
		double regionRho,regionP,regionT,regionK,regionOmega,regionMu_t,regionZ,regionZvar;
		// Assign specified values
		regionP=region.get_double("p");
		if (!FLAMELET) {
			if (region.get_double("rho").is_found) { // Rho is specified in the input file
				//  Check if Temperature is also specifid
				if (region.get_double("T").is_found) {
					cerr << "Both rho and T can't be specified in initial condition IC_" << ic+1 << endl;
					exit(1);
				}
				regionRho=region.get_double("rho");
				regionT=eos.T(regionP,regionRho);
			} else if (region.get_double("T").is_found) {
				regionT=region.get_double("T");
				regionRho=eos.rho(regionP,regionT);	
			} else {
				cerr << "Need to specify either rho or T in initial condition IC_" << ic+1 << endl;
				exit(1);
			}
		}
		regionK=region.get_double("k");
		regionOmega=region.get_double("omega");
		if (FLAMELET) {
			regionZ=region.get_double("Z");
			regionZvar=region.get_double("Zvar");
			double Chi=2.*regionOmega*regionZvar*rans.kepsilon.beta_star;
			regionRho=flamelet.table.get_rho(regionZ,regionZvar,Chi);
			regionT=flamelet.table.get_T(regionZ,regionZvar,Chi,false);
		}
		regionMu_t=regionRho*regionK/regionOmega;
		// If region is specified with a box method
		if (region.get_string("region")=="box") {
			// Loop the cells
			for (unsigned int c=0;c<grid.cellCount;++c) {
				// Check if the cell centroid is inside the box region
				if (withinBox(grid.cell[c].centroid,region.get_Vec3D("corner_1"),region.get_Vec3D("corner_2"))) {
					// Assign specified values
					grid.cell[c].p=regionP;
					grid.cell[c].T=regionT;
					grid.cell[c].rho=regionRho;
					grid.cell[c].v=regionV;
					if (TURBULENCE_MODEL!=NONE) {
						rans.cell[c].k=regionK;
						rans.cell[c].omega=regionOmega;
						rans.cell[c].mu_t=regionMu_t;
					}
					if (FLAMELET) {
						flamelet.cell[c].Z=regionZ;
						flamelet.cell[c].Zvar=regionZvar;
					}
				}
			}
		} else if (region.get_string("region")=="cylinder") {
			// Loop the cells
			for (unsigned int c=0;c<grid.cellCount;++c) {
				// Check if the cell centroid is inside the cylinder region
				Vec3D axisDirection=region.get_Vec3D("axisDirection");
				axisDirection=axisDirection.norm();
				if (withinCylinder(grid.cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"),axisDirection,region.get_double("height"))) {
					grid.cell[c].p=regionP;
					grid.cell[c].T=regionT;
					grid.cell[c].rho=regionRho;
					if (TURBULENCE_MODEL!=NONE) {
						rans.cell[c].k=regionK;
						rans.cell[c].omega=regionOmega;
						rans.cell[c].mu_t=regionMu_t;
					}
					if (FLAMELET) {
						flamelet.cell[c].Z=regionZ;
						flamelet.cell[c].Zvar=regionZvar;
					}
					// first component of the specified velocity is interpreted as the axial velocity
					// second component of the specified velocity is interpreted as the radial velocity
					// third component of the specified velocity is interpreted as the circumferential velocity
					Vec3D radialDirection=grid.cell[c].centroid-region.get_Vec3D("center");
					radialDirection=(radialDirection-radialDirection.dot(axisDirection)*axisDirection).norm();	
					grid.cell[c].v=		
						regionV[0]*axisDirection /* axial component */
						+regionV[1]*radialDirection /* radial component */
						+regionV[2]*radialDirection.cross(axisDirection); /* circumferential component */
				}
			}
		} else if (region.get_string("region")=="sphere") {
			// Loop the cells
			for (unsigned int c=0;c<grid.cellCount;++c) {
				// Check if the cell centroid is inside the sphere region
				if (withinSphere(grid.cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"))) {
					grid.cell[c].p=regionP;
					grid.cell[c].T=regionT;
					grid.cell[c].rho=regionRho;
					if (TURBULENCE_MODEL!=NONE) {
						rans.cell[c].k=regionK;
						rans.cell[c].omega=regionOmega;
						rans.cell[c].mu_t=regionMu_t;
					}
					if (FLAMELET) {
						flamelet.cell[c].Z=regionZ;
						flamelet.cell[c].Zvar=regionZvar;
					}
					// first component of the specified velocity is interpreted as the radial velocity
					// second and third components of the specified velocity are ignored
					grid.cell[c].v=regionV[0]*(grid.cell[c].centroid-region.get_Vec3D("center"));
				}
			}
		}
	}

	for (unsigned int c=0;c<grid.cellCount;++c) {
		for (unsigned int i=0;i<5;++i) grid.cell[c].update[i]=0.;
		if (FLAMELET) {
			double Chi=2.*rans.cell[c].omega*flamelet.cell[c].Zvar*rans.kepsilon.beta_star;
			flamelet.cell[c].mu=flamelet.table.get_mu(flamelet.cell[c].Z,flamelet.cell[c].Zvar,Chi);
		}
	}	
	
	for (unsigned int g=0;g<grid.ghostCount;++g) {
		for (unsigned int i=0;i<5;++i) grid.ghost[g].update[i]=0.;
	}
	
	if (TURBULENCE_MODEL!=NONE) for (unsigned int f=0;f<grid.faceCount;++f) rans.face[f].mu_t=0.;
	
	return;
}
