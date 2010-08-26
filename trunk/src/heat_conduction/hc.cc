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
#include "hc.h"

HeatConduction::HeatConduction (void) {
	// Empty constructor
	return;
}

void HeatConduction::initialize (void) {
	nVars=1;
	// Relative tolerance
	rtol=input.section("grid",gid).subsection("heatconduction").get_double("relativetolerance");
	// Absolute tolerance
	abstol=input.section("grid",gid).subsection("heatconduction").get_double("absolutetolerance");
	// Max linear solver iterations
	maxits=input.section("grid",gid).subsection("heatconduction").get_int("maximumiterations");

	mpi_init();
	material.set(gid);
	create_vars();
	apply_initial_conditions();
	set_bcs();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();
	petsc_init();

	return;
}

void HeatConduction::solve (int ts) {
	
	timeStep=ts;
	initialize_linear_system();
	assemble_linear_system();
	int nIter;
	double rNorm;
	petsc_solve(nIter,rNorm);
	if (Rank==0) cout << "\t" << nIter;
	update_variables();
	mpi_update_ghost_primitives();
	calc_cell_grads();
	mpi_update_ghost_gradients();
	
	return;
}

void HeatConduction::create_vars (void) {
	// Allocate variables
	// Default option is to store on cell centers and ghosts only
	T.allocate(gid);
	gradT.allocate(gid);
	update.allocate(gid);
	// qdot is only stored on bc faces
	qdot.cellStore=false; qdot.ghostStore=false; qdot.allocate(gid);
	return;
}

void HeatConduction::apply_initial_conditions (void) {
	
	double regionT;
	// Loop through each initial condition region and apply sequentially
	int count=input.section("grid",gid).subsection("IC",0).count;
	for (int ic=0;ic<count;++ic) {
		// Store the reference to current IC region
		Subsection &region=input.section("grid",gid).subsection("IC",ic);
		if (region.get_double("T").is_found) { // if T is specified in the input file
			regionT=region.get_double("T");
		} else {
			cerr << "T needs to be specified in initial condition IC_" << ic+1 << " of grid " << gid+1 << endl;
			exit(1);
		}
		
		// If region is specified with a box method
		if (region.get_string("region")=="box") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the box region
				if (withinBox(grid[gid].cell[c].centroid,region.get_Vec3D("corner_1"),region.get_Vec3D("corner_2"))) {
					// Assign specified values
					T.cell(c)=regionT;
				}
			}
		} else if (region.get_string("region")=="cylinder") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the cylinder region
				Vec3D axisDirection=region.get_Vec3D("axisdirection");
				axisDirection=axisDirection.norm();
				if (withinCylinder(grid[gid].cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"),axisDirection,region.get_double("height"))) {
					T.cell(c)=regionT;
				}
			}
		} else if (region.get_string("region")=="sphere") {
			// Loop the cells
			for (int c=0;c<grid[gid].cellCount;++c) {
				// Check if the cell centroid is inside the sphere region
				if (withinSphere(grid[gid].cell[c].centroid,region.get_Vec3D("center"),region.get_double("radius"))) {
					T.cell(c)=regionT;
				}
			}
		}
	}

	return;
}

void HeatConduction::calc_cell_grads (void) {
	vector<Vec3D> grad (3,0.);
	for (int c=0;c<grid[gid].cellCount;++c) {
		gradT.cell(c)=T.cell_gradient(c);
	 }
	return;
}

void HeatConduction::update_variables(void) {
	
	double residual=0.;
	double totalResidual;

	for (int c=0;c<grid[gid].cellCount;++c) {
		T.cell(c)+=update.cell(c);
		residual+=update.cell(c)*update.cell(c);
		
	} // cell loop
	
	MPI_Allreduce(&residual,&totalResidual,1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	if (timeStep==1) first_residual=sqrt(totalResidual);

	if (Rank==0) cout << "\t" << sqrt(totalResidual)/first_residual;
	
	return;
}


